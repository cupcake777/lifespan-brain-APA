# EasyColoc Tutorial

This tutorial is for users who want to run EasyColoc on real GWAS and QTL data. Start with the demo if this is your first time using the repository.

## 1. Install And Check The Environment

EasyColoc expects these command-line tools to be available: `Rscript`, `plink`, `bgzip`, and `tabix`.

Create the environment:

```bash
micromamba create -f environment.yml
micromamba activate easycoloc
```

Check the installation:

```bash
./easycoloc doctor
```

Run the repository smoke tests when you want a fuller check:

```bash
./easycoloc smoke
```

Run the small example check if you want a faster check of the CLI, report web server, and plot code:

```bash
bash examples/check_example.sh
```

It writes one plot to `examples/output/synthetic_locus_demo.pdf`. The input files in `examples/` are only for testing, not for real analysis.

## 2. Run The Demo First

The demo is the fastest way to see what a successful run looks like:

```bash
./easycoloc bootstrap-refs --demo ./demo_quickstart --run
```

Useful demo outputs:

| File or directory | What it is |
| --- | --- |
| `demo_quickstart/results/coloc_report.html` | static HTML summary |
| `demo_quickstart/results/all_colocalization_results.csv` | full coloc result table |
| `demo_quickstart/results/significant_colocalizations_PP4_*.csv` | thresholded coloc hits |
| `demo_quickstart/results/plots/` | locus plots |
| `demo_quickstart/results/rds/` | saved R objects for reuse |

If the demo runs, the software stack is basically working.

## 3. Create A Project For Your Data

Create a separate analysis folder instead of editing the repository defaults directly:

```bash
./easycoloc init /path/to/my_easycoloc_project
```

The new project contains:

| Path | Purpose |
| --- | --- |
| `config/global.yml` | output directory, reference files, thresholds, plot settings |
| `config/gwas.yml` | GWAS datasets, column names, genome build, sample size |
| `config/qtl.yml` | QTL metadata table, tabix file columns, QTL build |
| `data/` | a place to keep project input files |
| `results/` | recommended output directory |

## 4. Prepare The Inputs

For a real run, you need three groups of files:

| Input type | Examples |
| --- | --- |
| GWAS summary statistics | one or more GWAS files with SNP, chromosome, position, alleles, effect, SE, P value, and sample size information |
| QTL resources | tabix-indexed QTL `allPairs` files, optional `sigPairs` files, and a metadata table describing them |
| References | PLINK LD panel, allele frequency, reference genome, dbSNP, gene annotation, and optional recombination map |

The most important rule is that coordinates must be explicit. Record whether each GWAS and QTL file is `hg19` or `hg38`, and make the LD panel match the analysis build.

## 5. Optional: Build 1000 Genomes LD References

If you do not already have PLINK LD panels, EasyColoc can prepare them from 1000 Genomes data.

For `hg19`:

```bash
./easycoloc bootstrap-refs \
  --setup-1kg ./refs/1kg_phase3_hg19 \
  --build hg19 \
  --pop EAS \
  --chromosomes 1-22
```

For `hg38`:

```bash
./easycoloc bootstrap-refs \
  --setup-1kg ./refs/1kg_phase3_hg38 \
  --build hg38 \
  --pop EAS \
  --chromosomes 1-22
```

Use the population code that matches your GWAS as closely as possible.

## 6. Optional: Build GTEx Metadata

If your QTL data are local GTEx files, this helper creates a starting QTL config:

```bash
./easycoloc bootstrap-refs \
  --fetch-gtex-meta ./refs/gtex_meta \
  --gtex-eqtl-dir /path/to/gtex/eqtl \
  --gtex-sqtl-dir /path/to/gtex/sqtl
```

The command writes GTEx sample metadata, summary CSV files, and `qtl_gtex_generated.yml`.

## 7. Edit The Configs

Edit the config files in your project:

```text
/path/to/my_easycoloc_project/config/global.yml
/path/to/my_easycoloc_project/config/gwas.yml
/path/to/my_easycoloc_project/config/qtl.yml
```

Use `global.yml` for shared settings and reference paths. Use `gwas.yml` for GWAS file paths and column mapping. Use `qtl.yml` for QTL file metadata and tabix columns.

Before running, ask EasyColoc to show what files it expects:

```bash
./easycoloc refs \
  --global /path/to/my_easycoloc_project/config/global.yml \
  --gwas /path/to/my_easycoloc_project/config/gwas.yml \
  --qtl /path/to/my_easycoloc_project/config/qtl.yml \
  --include-qtl-files
```

Then validate the configuration:

```bash
./easycoloc doctor \
  --global /path/to/my_easycoloc_project/config/global.yml \
  --gwas /path/to/my_easycoloc_project/config/gwas.yml \
  --qtl /path/to/my_easycoloc_project/config/qtl.yml
```

Do not start a long run until `doctor` reports that the required files and tools are available.

## 8. Run The Pipeline

Use managed mode for normal work. It writes logs and runtime state files that make monitoring and debugging easier.

To run one configured GWAS x QTL pair without editing reduced config copies, list
the IDs and use `quick`:

```bash
./easycoloc quick /path/to/my_easycoloc_project --list
./easycoloc quick /path/to/my_easycoloc_project \
  --gwas-id MY_GWAS \
  --qtl-id MY_QTL
```

The command runs `doctor` for the selected pair first, keeps the normal GWAS
harmonization and reference checks, and writes to
`results/quick/MY_GWAS__MY_QTL/` by default. Use `--output-dir PATH` to override
that location.

For all configured datasets, run:

```bash
./easycoloc run --managed \
  --global /path/to/my_easycoloc_project/config/global.yml \
  --gwas /path/to/my_easycoloc_project/config/gwas.yml \
  --qtl /path/to/my_easycoloc_project/config/qtl.yml \
  --output-dir /path/to/my_easycoloc_project/results
```

## 9. Monitor Or Check A Run

During or after a run:

```bash
./easycoloc status /path/to/my_easycoloc_project/results
./easycoloc monitor /path/to/my_easycoloc_project/results
./easycoloc check /path/to/my_easycoloc_project/results
```

To save repeated monitor snapshots to a log file:

```bash
./easycoloc watch /path/to/my_easycoloc_project/results 60 logs/monitor/my_run.log
```

## 10. Review Results

Open the static report:

```text
/path/to/my_easycoloc_project/results/coloc_report.html
```

Or launch the local interactive report:

```bash
./easycoloc report-web /path/to/my_easycoloc_project/results
```

The main result files are:

| File | What to look for |
| --- | --- |
| `all_colocalization_results.csv` | all coloc tests and posterior probabilities |
| `significant_colocalizations_PP4_*.csv` | strongest candidate colocalizations |
| `all_susie_results.csv` | SuSiE fine-mapping summaries, if enabled |
| `plots/` | regional views of GWAS and QTL signals |
| `output_manifest.tsv` | list of generated output files |

## 11. Regenerate Plots

If the run already produced RDS bundles and you only changed plot settings, regenerate plots without rerunning the full coloc analysis:

```bash
Rscript tools/rerun_plots.R
```

## 12. Common Problems

| Symptom | What to check |
| --- | --- |
| `plink` is missing | Make sure the environment is active and `plink` is on `$PATH`. |
| `doctor` reports missing files | Run `./easycoloc refs --include-qtl-files` with the same config paths and fix the listed paths. |
| QTL files exist but no loci are tested | Check `qtl_info.file`, `allPairsTabixFilename`, `sigPairsTabixFilename`, and whether `sigPairs` filtering is too strict. |
| Results look shifted or empty | Recheck `hg19` versus `hg38` for GWAS, QTL, LD panel, allele frequency, annotation, and recombination map. |
| Interactive report does not open | Use the static `coloc_report.html`, or rerun `./easycoloc report-web RESULTS_DIR --no-open` and open the printed URL manually. |
