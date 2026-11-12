# --- paths ---
# Copy to paths.local.sh, set the values, then: set -a; source config/paths.local.sh; set +a
# APA_ROOT has the same layout on the cluster and locally:
#   raw/  ref/  data/<module>/  figures/  logs/
export APA_ROOT=/path/to/apa-pipeline
# Cluster tools and conda envs (SMR, PLINK2, MAGMA, LDSC)
export APA_HPC_HOME=/path/to/analysis_home
# Cohort inputs on the cluster (WGS/VCF, sample sheets; controlled access)
export APA_GROUP_HOME=/path/to/group_home
export APA_GROUP_SHARE=/path/to/group_share
export APA_GROUP_SHARE_ALT=/path/to/group_share_alt
export APA_PROTEOMICS_HOME=/path/to/proteomics_home
# Sentieon (00_geno joint calling)
export SENTIEON_HOME=/path/to/sentieon
export SENTIEON_LICENSE=/path/to/sentieon/license.lic
# Arial TTF folder for figures (optional)
export APA_FONT_DIR=/path/to/fonts
