# Module 07 handoff tables: 1 Mb loci, coloc evidence, provenance README.
# Reads and writes data/07_coloc; mediation input from data/04_qtl.
import csv
import json
import os
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# --- config ---
APA = Path(os.environ["APA_ROOT"])
ROOT = APA / "data/07_coloc"
README = ROOT / "README_handoff4.tsv"
PROJECT = Path(f"{os.environ['APA_HPC_HOME']}/project/work/apa")
FURONG = f"{os.environ['APA_HPC_HOME']}/project/work/apa"
# data/07_coloc names that differ from the names used on the source host
ORIG_NAME = {
    "figure_handoff.md": "figure.md",
    "chb_stage_eqtl_coloc_context.tsv": "chb_stage_eqtl_coloc.tsv",
    "gtex_brain_eqtl_coloc_tissues.tsv": "gtex_brain_eqtl_coloc_all_tissues.tsv",
    "smr_main_fig_source.tsv": "smr_main_figure_source.tsv",
    "Figure5c_magma_gsa_sum.tsv": "Figure5c_magma_gsa_summary.tsv",
    "assemble_main_figures_handoff.py": "assemble_main_figures.py",
    "assemble_main_figures_handoff2.py": "assemble_main_figures_v2.py",
}
GENES = {"SCAF1", "NRBP2", "ATP2B1", "MTRFR", "OGFOD2", "ARL3"}
GTEX_TISSUES = (
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
)


# --- helpers ---
def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return float("nan")


# --- context classification ---
signals = read_tsv(ROOT / "apa_signal_eqtl_classification.tsv")
gene_classes = {row["gene"]: row["gene_class"] for row in read_tsv(ROOT / "apa_gene_eqtl_classification.tsv")}
pair_classes = {(row["gene"], row["GWAS_ID"]): row for row in read_tsv(ROOT / "apa_gene_gwas_eqtl_classification.tsv")}
assert len(signals) == 104 and len(gene_classes) == 72

# Index windows are +/-500 kb; adjacent windows merge when index positions differ by <=1 Mb.
sites = sorted({(s["GWAS"], int(s["index_chrom"]), int(s["index_position"])) for s in signals})
locus_by_site = {}
counts = defaultdict(int)
last_position = {}
for gwas, chrom, position in sites:
    key = (gwas, chrom)
    if key not in last_position or position - last_position[key] > 1_000_000:
        counts[key] += 1
    locus_by_site[(gwas, chrom, position)] = f"{gwas}_chr{chrom}_L{counts[key]:03d}"
    last_position[key] = position

for signal in signals:
    signal["merged_1mb_locus_id"] = locus_by_site[(signal["GWAS"], int(signal["index_chrom"]), int(signal["index_position"]))]
write_tsv(ROOT / "signal_lookup_1mb.tsv", signals, list(signals[0]))

by_pair = defaultdict(list)
for signal in signals:
    by_pair[(signal["gene"], signal["GWAS"])].append(signal)

gene_rows = []
for (gene, gwas), group in sorted(by_pair.items()):
    pair = pair_classes[(gene, gwas)]
    chb = [number(s["chb_eqtl_max_pp4"]) for s in group]
    gtex = [number(s["gtex_max_pp4"]) for s in group]
    finite_chb = [v for v in chb if 0 <= v <= 1]
    finite_gtex = [v for v in gtex if 0 <= v <= 1]
    gene_rows.append({
        "gene": gene,
        "GWAS": gwas,
        "chb_max_pp4": max(finite_chb) if finite_chb else "",
        "gtex_brain_max_pp4": max(finite_gtex) if finite_gtex else "",
        "all_eqtl_resources_max_pp4": pair["eqtl_resource_max_pp4"],
        "has_eqtl_data": bool(finite_chb or finite_gtex),
        "apa_only_gene_gwas": pair["class"] == "APA-only",
        "apa_only_gene_all_gwas": gene_classes[gene] == "APA-only",
        "merged_1mb_locus_id": ";".join(sorted({s["merged_1mb_locus_id"] for s in group})),
    })
write_tsv(ROOT / "gene_gwas_classification_1mb.tsv", gene_rows, list(gene_rows[0]))

data_genes = {row["gene"] for row in gene_rows if row["has_eqtl_data"]}
apa_only_loci = {s["merged_1mb_locus_id"] for s in signals if gene_classes[s["gene"]] == "APA-only"}
assert len(gene_classes) == 72 and len(data_genes) <= 72 and len(apa_only_loci) <= 27
summary = [
    ("apa_gwas_colocalized_genes", 72, len(gene_classes)),
    ("genes_with_eqtl_data", 70, len(data_genes)),
    ("apa_only_genes", 28, sum(value == "APA-only" for value in gene_classes.values())),
    ("apa_only_merged_1mb_loci", 24, len(apa_only_loci)),
]
write_tsv(ROOT / "summary_old_new.tsv", [{"metric": k, "old": old, "corrected": new} for k, old, new in summary], ["metric", "old", "corrected"])


# --- representative coloc evidence ---
chb = {(row["signal_id"], row["chb_stage"]): row for row in read_tsv(ROOT / "chb_stage_eqtl_coloc_context.tsv")}
crossstage_path = ROOT / "representative_chb_crossstage.tsv"
if crossstage_path.exists():
    chb.update({(row["signal_id"], row["chb_stage"]): row for row in read_tsv(crossstage_path)})
gtex = {(row["signal_id"], row["tissue"]): row for row in read_tsv(ROOT / "gtex_brain_eqtl_coloc_tissues.tsv")}
evidence = []
for signal in signals:
    if signal["gene"] not in GENES:
        continue
    for stage in ("prenatal", "postnatal"):
        row = chb.get((signal["signal_id"], stage), {})
        evidence.append({
            "gene": signal["gene"], "signal_id": signal["signal_id"],
            "GWAS": signal["GWAS"], "apa_stage": signal["stage"],
            "apa_event": signal["apa_event"], "index_variant": signal["index_locus"],
            "apa_gwas_pp4": signal["apa_pp4"], "resource": "CHB", "resource_stage_or_tissue": stage,
            "stage_relation": "interaction" if signal["stage"] == "interaction" else ("matched" if stage == signal["stage"] else "opposite"),
            "eqtl_gwas_pp4": row.get("pp4", ""), "n_snps": row.get("n_snps", ""),
            "status": row.get("status", "not_tested"),
        })
    for tissue in GTEX_TISSUES:
        row = gtex.get((signal["signal_id"], tissue), {})
        evidence.append({
            "gene": signal["gene"], "signal_id": signal["signal_id"],
            "GWAS": signal["GWAS"], "apa_stage": signal["stage"],
            "apa_event": signal["apa_event"], "index_variant": signal["index_locus"],
            "apa_gwas_pp4": signal["apa_pp4"], "resource": "GTEx_v11", "resource_stage_or_tissue": tissue,
            "stage_relation": "adult_brain", "eqtl_gwas_pp4": row.get("pp4", ""),
            "n_snps": row.get("n_snps", ""), "status": row.get("status", "missing_source_row"),
        })
assert len(evidence) == 15 * sum(s["gene"] in GENES for s in signals)
write_tsv(ROOT / "representative_locus_eqtl_pp4.tsv", evidence, list(evidence[0]))

mediation = [row for row in read_tsv(APA / "data/04_qtl/04_mediation_results_module06.tsv") if row["gene_name"] == "ARL3"]
assert len(mediation) == 2
write_tsv(ROOT / "ARL3_mediation_rows.tsv", mediation, list(mediation[0]))


# --- figure source tables ---
for panel in ("a", "b", "c"):
    source = next(ROOT.glob(f"Figure6{panel}_*.json"))
    points = json.loads(source.read_text())["points"]
    fields = sorted({key for row in points for key in row})
    write_tsv(ROOT / f"Figure6{panel}_points.tsv", points, fields)

old_new = [
    ("eGenes_prenatal", 2113, 4735, "permutation-significant eGenes"),
    ("eGenes_postnatal", 9713, 13929, "permutation-significant eGenes"),
    ("APA_eQTL_coloc_prenatal_PP4_ge_0.7", "165/666", "186/665", "APA-locus denominator"),
    ("APA_eQTL_coloc_postnatal_PP4_ge_0.7", "576/2231", "655/2229", "APA-locus denominator"),
    ("mediation_pairs_PP_H4_ge_0.8", 658, 772, "137+521 versus 164+608"),
    ("APA_only_genes", 28, summary[2][2], "gene across matched GWAS"),
    ("APA_only_merged_1mb_loci", 24, summary[3][2], "within GWAS and chromosome"),
]
write_tsv(ROOT / "old_new_comparison.tsv", [{"metric": a, "old": b, "corrected": c, "definition": d} for a, b, c, d in old_new], ["metric", "old", "corrected", "definition"])


# --- source manifest ---
remote_context = {
    "chb_stage_eqtl_coloc_context.tsv", "gtex_brain_eqtl_coloc_tissues.tsv",
    "apa_signal_eqtl_classification.tsv", "apa_gene_gwas_eqtl_classification.tsv",
    "apa_gene_eqtl_classification.tsv", "apa_only_summary.tsv", "analysis_summary.md",
}
module07_disease = {"07_internal_eqtl_integration_contract.tsv.gz", "internal_summary.tsv", "04_disease_chains.tsv.gz", "05_evidence_tiers.tsv.gz"}
generated = {"signal_lookup_1mb.tsv", "gene_gwas_classification_1mb.tsv", "summary_old_new.tsv", "representative_locus_eqtl_pp4.tsv", "ARL3_mediation_rows.tsv", "old_new_comparison.tsv"}
job = os.environ.get("SLURM_JOB_ID", "NA")
manifest = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or path == README:
        continue
    rel = path.relative_to(ROOT).as_posix()
    old = ORIG_NAME.get(rel, rel)
    host, source, producer, job_id = "sk-bio", str(PROJECT / "pub" / old), "source file", "NA"
    if rel in ("figure_handoff.md", "results.md"):
        source, producer = str(PROJECT / "pub" / old), "manuscript register"
    elif rel in remote_context:
        host, source, producer, job_id = "furong", f"{FURONG}/14_context/{old}", "14_context/apa_only_eqtl_coloc_v2.R", "46323208"
        if rel == "analysis_summary.md":
            producer, job_id = "14_context/analysis_summary.md", "NA"
    elif rel == "representative_chb_crossstage.tsv":
        host, source, producer, job_id = "furong", f"{FURONG}/14_context/representative_chb_crossstage.tsv", "14_context/representative_chb_crossstage.R", "46325185"
    elif rel == "eqtl_coloc.tsv":
        host, source, producer, job_id = "furong", f"{FURONG}/07_locus_colocalization_smr/func/results/eqtl_coloc.tsv", "func/scripts/02c_eqtl_coloc.R", "46324566"
    elif rel in module07_disease:
        host, source, producer, job_id = "furong", f"{FURONG}/05_func/{'work/coloc' if rel == 'internal_summary.tsv' else 'results'}/{rel}", "05_func/scripts/06_disease.py", "46325111"
    elif rel == "Figure5c_magma_gsa_sum.tsv":
        source, producer = str(PROJECT / "10_magma/enrichment/magma_gsa_summary.tsv"), "prior Figure 5C MAGMA source; replaced by disease gene matrix"
    elif rel.startswith("Figure5c_"):
        source, producer, job_id = str(path), "05_disease_gene_matrix.py", "84148"
    elif rel.startswith("Figure6") and rel.endswith(".json"):
        panel_source = {
            "a": "EAS_SCZ_interaction_chr12_89410263_rs12578201_ATP2B1_locuscompare.json",
            "b": "EAS_SCZ_postnatal_chr12_122799848_rs77908463_MTRFR_locuscompare.json",
            "c": "EAS_SCZ_postnatal_chr10_102760520_rs12571643_ARL3_locuscompare.json",
        }
        source, producer = str(PROJECT / "04_coloc/results/plot_data" / panel_source[path.name[7]]), "04_coloc/05_plot_locuscompare.R"
    elif rel.startswith("Figure6") and rel.endswith("_points.tsv"):
        source, producer, job_id = str(next(ROOT.glob(f"Figure6{path.name[7]}_*.json"))), "build_handoff.py", job
    elif rel.startswith("smr_"):
        source, producer = str(PROJECT / "06_smr/results/plots" / old), "06_smr/03_plot_smr_style.R"
    elif rel == "Figure6d_ARL3_smr.tsv":
        source, producer, job_id = str(PROJECT / "06_smr/results/smr_all_results.tsv"), "06_smr/03_plot_smr_style.R; handoff ARL3 subset", "84137"
    elif rel in generated:
        source, producer, job_id = "; ".join(f"{FURONG}/14_context/{name}" for name in ("apa_signal_eqtl_classification.tsv", "chb_stage_eqtl_coloc.tsv", "gtex_brain_eqtl_coloc_all_tissues.tsv")), "build_handoff.py", job
        if rel == "ARL3_mediation_rows.tsv":
            source = str(PROJECT / "11_eqtl/results/04_mediation_results.tsv")
        elif rel == "old_new_comparison.tsv":
            source = "; ".join((str(PROJECT / "pub/feedback.md"), str(PROJECT / "pub/.mem"), str(ROOT / "summary_old_new.tsv")))
        elif rel == "representative_locus_eqtl_pp4.tsv" and crossstage_path.exists():
            source += f"; {FURONG}/14_context/representative_chb_crossstage.tsv"
    elif rel.startswith(("Figure5_", "Figure6_provisional")):
        source, producer, job_id = "assemble_main_figures.py", "assemble_main_figures.py", "84137"
        if rel.startswith("Figure5_"):
            source, producer, job_id = "assemble_main_figures_v2.py", "assemble_main_figures_v2.py", "84149"
    elif rel == "assemble_main_figures_handoff.py":
        source, producer = "assemble_main_figures.py", "copy with handoff output path"
    elif rel == "assemble_main_figures_handoff2.py":
        source, producer = "assemble_main_figures_v2.py", "copy with handoff output path"
    manifest.append({"filename": rel, "source_path": source, "producer_script": producer,
                     "job_id": job_id, "modified_time": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds"),
                     "compute_host": host})
write_tsv(README, manifest, ["filename", "source_path", "producer_script", "job_id", "modified_time", "compute_host"])
