#!/usr/bin/env bash
# Copy to config/paths.local.sh and edit for the local environment.

export APA_CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export APA_ROOT="/path/to/apa-data"
export APA_FIGURE_DIR="${APA_CODE_ROOT}/outputs/figures"

# Optional: directory containing arial.ttf, arialbd.ttf, ariali.ttf and arialbi.ttf.
# R and Python plotting code fall back to an available sans-serif font when unset.
export APA_FONT_DIR=""

# Required only for cluster-scale analyses that call site-installed tools.
export APA_HPC_HOME="/path/to/cluster/software-and-work"
export APA_GROUP_HOME="/path/to/group/storage"
export APA_GROUP_SHARE="/path/to/group/share"
export APA_PROTEOMICS_HOME="/path/to/proteomics/source"
export ANNO_DIR="/path/to/external/annotation/resources"

# External analysis tools are installed separately rather than vendored here.
export EASYCOLOC_ROOT="/path/to/easycoloc"
