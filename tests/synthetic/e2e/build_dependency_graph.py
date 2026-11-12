#!/usr/bin/env python3
# Static dependency graph of the release scripts: file names read and written, matched by basename
# (glob aware), plus the consumed files that no release script produces.
# Writes dependency_graph.tsv and missing_links.tsv to $APA_ROOT/data/tests/.
# Heuristic: a token's role is the nearest I/O keyword on its line; run-time-built names are not resolved except through OVERRIDES.

import fnmatch
import glob
import os
import re

import pandas as pd

# --- helpers ---
def is_glob(b):
    return any(c in b for c in '*?[')


# --- paths ---
HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
APA_ROOT = os.environ['APA_ROOT']
OUT = os.path.join(APA_ROOT, 'data', 'tests')
os.makedirs(OUT, exist_ok=True)
RESULTS_MANIFEST = os.path.join(APA_ROOT, 'data', '05_gxs', 'RESULTS_MANIFEST.tsv')

# --- patterns and tables ---
DATA_RE = re.compile(r'''(?<![\w/.\-])((?:[\w.\-${}/*]+/)*[\w.\-${}*]+\.(?:tsv|txt|csv|rds|rda|RData|rdata|gz|bed|json|xlsx|xls|pdf|png|svg|md|obo|gmt|zip|pgen|pvar|psam|vcf|meta|fa|fai|parquet|pkl|h5|bim|fam|bgen|sumstats|annot|ld|hsq|log|yml|yaml|cfg|config|sample|final|matrix|bim|xz|bgz|tab|res|out|list|mat|npz|tbl|dat))(?![\w])''')
GOOD_EXT = ('tsv', 'txt', 'csv', 'rds', 'rda', 'RData', 'rdata', 'gz', 'bed', 'json', 'xlsx', 'xls', 'pdf', 'png',
            'svg', 'md', 'obo', 'gmt', 'zip', 'pgen', 'pvar', 'psam', 'vcf', 'meta', 'fa', 'fai', 'parquet', 'pkl',
            'h5', 'bim', 'fam', 'sumstats', 'annot', 'ld', 'hsq', 'yml', 'yaml', 'cfg', 'config', 'sample', 'final',
            'out', 'sbatch', 'slurm')
METHOD_PREFIX = re.compile(r'^(as|np|math|pd|os|plt|sns|json|sys|dt|tools|utils|ax|fig|df|self|obj|x|res|fit|mod|'
                           r'scipy|stats|glm|gam|base|data|mat|nrow|ldscore|ggplot|str|file|tmp|out|df2|dat|cat|sub|'
                           r'rbind|table)\.')
WRITE_KW = re.compile(r'(fwrite|write[._\w]*|save[A-Za-z_.]*|to_csv|to_excel|to_parquet|ggsave|pdf\s*\(|png\s*\(|'
                      r'svg\s*\(|jpeg\s*\(|savefig|dump|\btee\b|sink\s*\(|--out[\w-]*|\s-o\s|>\s*\S|>>|output|outfile|'
                      r'out_dir|outdir|(?-i:\bOUT\w*)|(?-i:\bDEST\w*)|export|gzip|bgzip|writeLines|create_|\bcp\s|\bmv\s|'
                      r'mkdir|file\s*=\s*out|\bcat\s*\(|>\s*$|write_table|write_tsv|json\.dump|\.open\(\s*["\']w|'
                      r'\bsprintf\s*\(|\bpaste0?\s*\(|\.to_csv|\.to_parquet)', re.I)
READ_KW = re.compile(r'(\bfread|read[._\w]*|readRDS|\bload\s*\(|source\s*\(|import\b|open\s*\(|readLines|scan\s*\(|'
                     r'gzfile|zcat|\bzgrep|(?<![<\-])<\s*(?![-<])\S|--in[\w-]*|--sumstats|--bfile|--pheno|--gwas|'
                     r'--input|--annot|--ref|--geno|--covariates|--bed|\bin\s*=|(?-i:\bIN\w*)|input|Sys\.getenv|'
                     r'os\.environ|file\.exists|list\.files|grep|awk|head\s|tail\s|tabix|bcftools|--file|\bcol_|'
                     r'(?-i:\bREF\w*)|\bread\b)', re.I)
WRITE_CALL = re.compile(r'(fwrite|write\.|write_|writeLines|saveRDS|save\(|ggsave|to_csv|to_excel|to_parquet|savefig|'
                        r'\.dump\(|pdf\(|png\(|svg\(|sink\(|\s>\s|--out|\btee\b|cp\s)')
WRITE_LINE = re.compile(r'(fwrite|write\w*|saveRDS|ggsave|to_csv|to_excel|to_parquet|savefig|dump\(|'
                        r'\.open\(\s*["\']w|sink\(|save\w*\()', re.I)
EXT_RES = re.compile(r'\.gmt$|^c[1-9]\.|GO_|Allen_Brain|Human_Phenotype|HGNC|MitoCarta|PolyA|POSTAR|RBPDB|GTEx|kegg|goslim|SynGO|syngo|BrainSpan|brainspan|gencode|hg19|hg38|\.chain|miRNA|clip\.csv|cibersortx|Dapars2_config|refseq', re.I)
SAMPLE_RE = re.compile(r'\.meta$|sample\.(Adult|Prenatal|Postnatal)|samples?\.(pre|post)natal|\.sample$|sample_list', re.I)
HPC_RE = re.compile(r'parquet$|susie|cis_qtl|signif_pairs|\.ldscore|\.l2\.|\.gsa\.out|sumstats|\.annot$|_summary\.(tsv|csv)|coloc_results|\.pairs|_pairs|\.bim$|\.pvar$|\.psam$|\.pgen$|\.fam$', re.I)
ORDER = {'00_geno': 0, '01_quant': 1, '02_traj': 2,
         '03_omics': 3, '04_qtl': 4, '05_gxs': 5,
         '06_disease': 6, '07_coloc': 7, 'figures': 8, 'tests': 9}
TRACK_KEYS = {'magma': 'MAGMA summary (coloc/MAGMA track)', 'coloc': 'coloc summary (coloc track)',
              'sldsc': 'S-LDSC summary (S-LDSC 2-stage BH track)', 'gating': 'S4 gating table (S4 track)',
              's22': 'S2 regulatory-element table (S2 track)', 'elements': 'S2 regulatory-element table (S2 track)'}
# Manual resolutions, each checked by reading the cited line(s). value = (category, note, writer or '')
OVERRIDES = {
    'significant_unique_trait_phenotype_PP4_0.7.csv': ('vendored_test_artifact', 'EasyColoc tools/checks smoke-test output', ''),
    'coloc_report_data.json': ('vendored_test_artifact', 'EasyColoc tools/checks smoke-test output', ''),
    'overview_input_audit.tsv': ('resolved_writer_in_release', 'writer: 06 TWAS 07_multi_overview.py (to_csv via OVERVIEW_AUDIT_PATH); parser miss', '06_disease/twas/07_multi_overview.py'),
    'interaction_qtl_bridge_input_audit.tsv': ('resolved_writer_in_release', 'writer: 06 TWAS 08_inter_bridge.py (INPUT_AUDIT_PATH); no downstream consumer', '06_disease/twas/08_inter_bridge.py'),
    'args.config': ('false_positive', 'argparse attribute, not a file', ''),
    'figure4_dynamic_input_screen.tsv': ('candidate_missing_producer', 'read by QTL/plot_figure4_dynamic.py; no writer in release', ''),
    'figure4_dynamic_input_wig_profile_qc.tsv': ('candidate_missing_producer', 'read by QTL/plot_figure4_dynamic.py; no writer in release', ''),
    'figure4_dynamic_input_wig_profile_samples.tsv': ('candidate_missing_producer', 'read by QTL/plot_figure4_dynamic.py; no writer in release', ''),
    'genotype.prenatal.pvar': ('resolved_templated_writer', 'pre_input/genotype.%s written by plink2 --make-pgen in 04 eqtl/01_preprocess.R (sprintf)', '04_qtl/eqtl/01_preprocess.R'),
    'genotype.postnatal.pvar': ('resolved_templated_writer', 'pre_input/genotype.%s written by plink2 --make-pgen in 04 eqtl/01_preprocess.R (sprintf)', '04_qtl/eqtl/01_preprocess.R'),
    'genotype.psam': ('candidate_missing_producer', 'read by 07 func 12_plot_hid1_mechanism.py; no writer of this prefix in release', ''),
    'enrichment_GO_shared.tsv': ('candidate_missing_producer|being_reconstructed_by_other_track', 'MAGMA GO summary; run_enrichment.R writes only prenatal_vs_postnatal', ''),
    'enrichment_GO_postnatal_only.tsv': ('candidate_missing_producer|being_reconstructed_by_other_track', 'MAGMA GO summary; run_enrichment.R writes only prenatal_vs_postnatal', ''),
    'hg38_autr_candidate_intervals.fa': ('candidate_missing_producer', 'derived FASTA read by 07 func 12_autr_mirna_rbp.py; no writer in release', ''),
    'Intron_ARE.csv': ('documented_external', '07 func reference set (config.yml)', ''),
    'EAS_interaction.tsv': ('candidate_missing_producer|being_reconstructed_by_other_track', 'S-LDSC interaction input named in 04 eqtl 04-3_sldsc_interaction.sh; no writer', ''),
    'magma_traj_summary.tsv': ('present_in_worktree_no_release_producer|being_reconstructed_by_other_track', 'MAGMA trajectory summary; server run per registry; hand-off copy only', ''),
    'cibersortx_mixture_prenatal.tsv': ('external_tool_output_not_run', '07 func 03_deconvolution.R records cibersortx_status not_run', ''),
    'cibersortx_mixture_postnatal.tsv': ('external_tool_output_not_run', '07 func 03_deconvolution.R records cibersortx_status not_run', ''),
    'samples.prenatal.txt': ('withheld_sample_list', 'per-sample list; WITHHELD apa.sample.*.final category', ''),
    'samples.postnatal.txt': ('withheld_sample_list', 'per-sample list; WITHHELD apa.sample.*.final category', ''),
    'sample.Adult.final': ('withheld_sample_list', 'per-sample list; WITHHELD apa.sample.*.final category', ''),
    'annotations.xlsx': ('documented_external', 'SynGO annotation workbook (external; cited in 02 followup)', ''),
    'genes.xlsx': ('documented_external', 'reference workbook (external; not in release)', ''),
    'ontologies.xlsx': ('documented_external', 'reference workbook (external; not in release)', ''),
    'refseq_to_gene.txt': ('documented_external', 'DaPars2 reference table (external; not in release)', ''),
    'chrname.txt': ('documented_external', 'chromosome-name reference (external)', ''),
    'Chrname.txt': ('documented_external', 'chromosome-name reference (external)', ''),
    'Brain_Frontal_Cortex_BA9.v11.covariates.txt': ('documented_external', 'GTEx v11 covariates (EXTERNAL_INPUTS)', ''),
    '1-s2.0-S0896627318302368-mmc5.xls': ('documented_external', 'published supplement (EXTERNAL_INPUTS)', ''),
    '_snp151_hash_table.json': ('documented_external', 'dbSNP hash table (EasyColoc tool)', ''),
    'gencode.v40.tx2gene.tsv': ('candidate_missing_producer', 'runtime-confirmed in dry run: 02 pipeline/05_func.R stops with this file missing; derived from external GENCODE v40 GTF; no producer in release', ''),
    'read.csv': ('false_positive', 'R function name, not a file', ''),
    'MANIFEST.tsv': ('release_metadata', 'self-check target of tests/check_portable_tree.sh', ''),
    'environment.yml': ('release_metadata', 'self-check target of tests/check_portable_tree.sh', ''),
    'qc1.fam': ('tool_output_plink', 'PLINK --out stem produced by 04_preQC.sh', ''),
    'qc1_clean.fam': ('tool_output_plink', 'PLINK --out stem produced by 04_preQC.sh', ''),
    'qc3.fam': ('tool_output_plink', 'PLINK --out stem produced by 04_preQC.sh', ''),
    '-Dapars2_config.txt': ('shell_token_artifact', 'leading-dash token; DaPars2 config from generate_configure_for_Dapars2.py', ''),
    'DEMO.md': ('vendored_test_artifact', 'EasyColoc demo fixture', ''),
    'input.tsv': ('vendored_test_artifact', 'EasyColoc test fixture', ''),
}

# --- scan release scripts for file-name tokens ---
scripts = []
for root, _, files in os.walk(REL):
    if '/.git' in root or '__pycache__' in root:
        continue
    for f in files:
        rel = os.path.relpath(os.path.join(root, f), REL)
        if rel.startswith('tests' + os.sep + 'synthetic'):
            continue  # e2e test tooling is not part of the analysis graph
        if rel.endswith(('.R', '.r', '.py', '.sh', '.slurm', '.sbatch', '.bash')):
            scripts.append(rel)
scripts.sort()
text = {s: open(os.path.join(REL, s), encoding='utf-8', errors='replace').read() for s in scripts}
has_write = {s: bool(WRITE_CALL.search(t)) for s, t in text.items()}
rows = []
for s in scripts:
    for i, line in enumerate(text[s].split('\n'), 1):
        if line.lstrip().startswith('#'):
            continue
        for m in DATA_RE.finditer(line):
            tok = m.group(1)
            b = tok.split('/')[-1]
            if b.endswith(('.R', '.r', '.py', '.sh')):
                continue
            ext = b.rsplit('.', 1)[-1]
            if ext not in GOOD_EXT and not b.endswith(('.tsv.gz', '.txt.gz', '.bed.gz', '.gz')):
                continue
            if '/' not in tok and METHOD_PREFIX.match(b):
                continue
            # role from the nearest I/O keyword in the 80 characters before the token
            win = line[max(0, m.start() - 80):m.start()]
            best = None
            for kind, rx in (('write', WRITE_KW), ('read', READ_KW)):
                for km in rx.finditer(win):
                    if best is None or km.end() >= best[2]:
                        best = (kind, km.group(0), km.end())
            role = best[0] if best else None
            # no keyword: assignment target name decides
            if role is None:
                am = re.match(r'\s*([A-Za-z_][\w.]*)\s*(?:=|<-|:=)', line)
                if am:
                    if re.search(r'(?i)(^|_)(out|output|save|dest|write|outdir|outfile)(_|$)|^OUT', am.group(1)):
                        role = 'write'
                    else:
                        role = 'read'
            if role is None:
                w = bool(WRITE_KW.search(line))
                r_ = bool(READ_KW.search(line))
                role = 'write' if (w and not r_) else ('read' if (r_ and not w) else 'unknown')
            if role == 'unknown':
                erole = 'write?' if has_write[s] else 'read?'
            else:
                erole = role
            rows.append(dict(script=s, line=i, erole=erole, basename=b, token=tok))
E = pd.DataFrame(rows)

# --- variable-path writes: path literal bound to a variable, passed to a write call in the same file ---
extra, drop = [], set()
for s_ in scripts:
    lines = text[s_].split('\n')
    bound = {}
    for i, line in enumerate(lines, 1):
        if line.lstrip().startswith('#'):
            continue
        m = re.match(r'\s*([A-Za-z_][\w.]*)\s*(?:=|<-|:=)', line)
        if not m:
            continue
        for mm in DATA_RE.finditer(line):
            b = mm.group(1).split('/')[-1]
            if b and '{' not in b and '$' not in b and not b.endswith(('.R', '.r', '.py', '.sh')):
                bound.setdefault(m.group(1), set()).add(b)
    if not bound:
        continue
    for i, line in enumerate(lines, 1):
        if line.lstrip().startswith('#') or not WRITE_LINE.search(line):
            continue
        for var, bases in bound.items():
            if re.search(r'(?<![\w.])' + re.escape(var) + r'(?![\w])', line):
                for b in bases:
                    extra.append(dict(script=s_, line=i, erole='write', basename=b, token=b))
                    drop.add((s_, b))
if extra:
    E = pd.concat([E, pd.DataFrame(extra)], ignore_index=True)
    E = E[~E.apply(lambda r: r.erole.startswith('read') and (r.script, r.basename) in drop, axis=1)]
    E = E.drop_duplicates(subset=['script', 'basename', 'erole'])
prod = E[E.erole.isin(['write', 'write?'])]
prodmap = prod.groupby('basename')['script'].apply(lambda x: sorted(set(x)))
prod_keys = set(prodmap.index)
Rc = E[E.erole.isin(['read', 'read?'])].copy()

# --- external and withheld documentation text ---
ext_text = ''
for f in glob.glob(os.path.join(REL, '_manifests', '*EXTERNAL_INPUTS.tsv')) + \
        glob.glob(os.path.join(REL, '_manifests', '*WITHHELD.tsv')):
    ext_text += ' ' + open(f, encoding='utf-8', errors='replace').read()

# --- working-tree presence (read-only) ---
wt = {}
for root, _, files in os.walk(APA_ROOT):
    if root.startswith(REL) or '/.git' in root or '/.claude' in root:
        continue
    for f in files:
        wt.setdefault(f, []).append(os.path.relpath(os.path.join(root, f), APA_ROOT))

# --- dependency_graph.tsv: match each read token to a producer ---
links = []
unmatched = []
for _, r in Rc.iterrows():
    b = r.basename
    # producer lookup: exact name, then glob either way
    prods, how = [], ''
    if b in prod_keys:
        prods, how = prodmap[b], 'exact'
    else:
        hits = [k for k in prod_keys if fnmatch.fnmatch(k, b)] if is_glob(b) else []
        if hits:
            prods, how = sorted(set(sum([prodmap[k] for k in hits], []))), 'glob'
        else:
            for k in prod_keys:
                if is_glob(k) and fnmatch.fnmatch(b, k):
                    prods, how = prodmap[k], 'glob'
                    break
    prods = [p for p in prods if p != r.script]
    cons_mod = r.script.split('/')[0]
    # documented external or withheld reference: name appears in the EXTERNAL/WITHHELD manifests
    stem = re.sub(r'(\.gz|\.tsv|\.txt|\.csv|\.rds|\.bed|\.json|\.xlsx|\.pdf|\.png|\.svg|\.md|\.zip)+$', '', b)
    doc_ext = len(stem) >= 6 and '*' not in b and not stem.startswith('$') and ((b in ext_text) or (stem in ext_text))
    if prods:
        for p in prods:
            links.append(dict(producer=p, producer_module=p.split('/')[0], file=r.basename, consumer=r.script,
                              consumer_module=cons_mod, match=how,
                              backward_order=ORDER.get(p.split('/')[0], -1) > ORDER.get(cons_mod, -1),
                              link_status='linked'))
    elif doc_ext:
        links.append(dict(producer='(external/withheld, documented)', producer_module='external', file=r.basename,
                          consumer=r.script, consumer_module=cons_mod, match='documented_external',
                          backward_order=False, link_status='external_documented'))
    else:
        unmatched.append(dict(producer='NONE', producer_module='', file=r.basename, consumer=r.script,
                              consumer_module=cons_mod, match='unmatched', backward_order=False,
                              link_status='unmatched'))
DG = pd.DataFrame(links + unmatched).drop_duplicates()
DG = DG.sort_values(['link_status', 'producer', 'file', 'consumer']).reset_index(drop=True)
DG.to_csv(os.path.join(OUT, 'dependency_graph.tsv'), sep='\t', index=False)

# --- missing_links.tsv: unmatched consumed files, categorised ---
miss = pd.DataFrame(unmatched)
cons_files = set(Rc.basename)
mrows = []
if len(miss):
    for b, g in miss.groupby('file'):
        cons = sorted(set(g.consumer))
        wt_paths = wt.get(b, [])
        if not wt_paths and is_glob(b):
            wt_paths = [p for k in wt for p in wt[k] if fnmatch.fnmatch(k, b)][:2]
        cat, note, writer = 'candidate_missing_producer', '', ''
        other = [v for k, v in TRACK_KEYS.items() if k in b.lower()]
        if b in OVERRIDES:
            cat, note, writer = OVERRIDES[b]
        elif '{' in b or '$' in b or '%' in b:
            cat, note = 'templated_name', 'name built from a variable or template; producer depends on that variable'
        elif any(c.startswith('figures/') for c in cons) and not wt_paths:
            cat, note = 'figure_data_producer_absent', 'figure-data table consumed by figures/*.R; producer not in release'
        elif EXT_RES.search(b):
            cat, note = 'documented_external_reference_pattern', 'reference resource matched by name pattern (not shipped; see EXTERNAL_INPUTS)'
        elif SAMPLE_RE.search(b):
            cat, note = 'individual_level_input_not_shipped', 'sample metadata or sample list (WITHHELD category; not shipped)'
        elif HPC_RE.search(b) and not wt_paths:
            cat, note = 'hpc_tool_output_not_in_release', 'tensorQTL/SuSiE/LDSC/MAGMA/coloc-type output; server-side producer not in release'
        elif any(c.startswith('figures/') for c in cons):
            cat = 'figure_data_in_worktree_no_producer' if wt_paths else 'figure_data_producer_absent'
            note = ('exists only in working-tree hand-off ' + wt_paths[0] + '; no release producer') if wt_paths \
                else 'figure-data table consumed by figures/*.R; producer not in release'
        elif wt_paths:
            cat, note = 'present_in_worktree_no_release_producer', 'exists only in working-tree hand-off ' + wt_paths[0]
        if other and cat.startswith(('candidate', 'present', 'figure')) and 'being_reconstructed' not in cat:
            cat = cat + '|being_reconstructed_by_other_track'
            note = (note + '; ' if note else '') + other[0]
        mrows.append(dict(link_type='missing_producer', file=b, final_category=cat, consumer_scripts='; '.join(cons)[:600],
                          consumer_modules=','.join(sorted(set(c.split('/')[0] for c in cons))), n_consumers=len(cons),
                          working_tree_location=(wt_paths[0] if wt_paths else ''), writer_in_release=writer, notes=note))

# --- orphan outputs: strong writes, no static consumer, script not mapped to a manuscript item ---
man_scripts = set()
md = pd.read_csv(os.path.join(REL, 'MANIFEST.tsv'), sep='\t', dtype=str).fillna('')
for _, r in md.iterrows():
    if r.manuscript_items.strip() and r.status != 'superseded_kept_for_provenance':
        man_scripts.add(r.release_path.replace('release/', '', 1))
rm_names = set()
rmd = pd.read_csv(RESULTS_MANIFEST, sep='\t', dtype=str).fillna('')
for v in rmd.produced_by.str.replace('release/', '', regex=False):
    for m in re.finditer(r'([\w/\.\-]+\.(?:R|r|py|sh|slurm|sbatch))', v):
        rm_names.add(os.path.basename(m.group(1)))
strong = E[(E.erole == 'write') & ~E.script.str.startswith(('figures/', 'tests/'))][['script', 'basename']].drop_duplicates()
strong = strong[(strong.basename != '') & ~strong.basename.isin(cons_files)]
orph = strong.groupby('script').basename.apply(lambda x: '; '.join(sorted(set(x)))[:300]).reset_index()
for _, r in orph.iterrows():
    if r.script in man_scripts or os.path.basename(r.script) in rm_names or r.script.startswith('figures/'):
        continue
    mrows.append(dict(link_type='orphan_output_script', file=r.basename, final_category='no_static_consumer_and_no_manuscript_map',
                      consumer_scripts='', consumer_modules='', n_consumers=0, working_tree_location='',
                      writer_in_release='', notes='script ' + r.script + '; unconsumed outputs listed in file column'))
ML = pd.DataFrame(mrows)
ML.to_csv(os.path.join(OUT, 'missing_links.tsv'), sep='\t', index=False)
print('scripts', len(scripts), 'read/write tokens', len(E), 'links', len(links), 'unmatched rows', len(unmatched),
      'missing names', (ML.link_type == 'missing_producer').sum(), 'orphan rows', (ML.link_type == 'orphan_output_script').sum())
print(ML[ML.link_type == 'missing_producer'].final_category.value_counts().to_string())
