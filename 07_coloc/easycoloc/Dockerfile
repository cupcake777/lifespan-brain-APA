# EasyColoc Docker image
#
# Build:
#   docker build -t easycoloc:latest .
#
# Run:
#   docker run -v /path/to/project:/work -v /path/to/ref:/ref:ro easycoloc:latest run --managed \
#     --global /work/config/global.yml --gwas /work/config/gwas.yml --qtl /work/config/qtl.yml --output-dir /work/results
#
# Interactive:
#   docker run -it --entrypoint bash -v /path/to/project:/work easycoloc:latest

FROM rocker/verse:4.3.2

LABEL maintainer="EasyColoc Team"
LABEL description="Fast and easy-to-use colocalization analysis pipeline"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# System Dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools
    build-essential \
    cmake \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Genomics tools
    plink \
    tabix \
    bgzip \
    bedtools \
    samtools \
    bcftools \
    # Python
    python3 \
    python3-pip \
    python3-venv \
    # Utilities
    wget \
    curl \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# R Package Dependencies
# Install Bioconductor manager first
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'data.table', \
    'dplyr', \
    'yaml', \
    'parallel', \
    'glue', \
    'coloc', \
    'susieR', \
    'ggplot2', \
    'ggpubr', \
    'vroom', \
    'forcats', \
    'purrr', \
    'ggrepel', \
    'jsonlite', \
    'tidyr', \
    'magrittr', \
    'tibble', \
    'readr', \
    'stringr' \
    ), repos='https://cloud.r-project.org', Ncpus=4)"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'ensembldb', \
    'GenomicRanges', \
    'GenomicFeatures', \
    'AnnotationDbi', \
    'rtracklayer', \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'TxDb.Hsapiens.UCSC.hg38.knownGene', \
    'VariantAnnotation', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'BSgenome.Hsapiens.UCSC.hg19' \
    ), Ncpus=4, ask=FALSE, update=FALSE)"

# =============================================================================
# Python Dependencies
# =============================================================================
# Create virtual environment
RUN python3 -m venv /opt/easycoloc_py
ENV PATH="/opt/easycoloc_py/bin:$PATH"

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install --no-cache-dir \
    adjustText \
    biopython \
    gtfparse \
    liftover \
    matplotlib \
    numpy \
    pandas \
    pyensembl \
    pysam \
    scikit-allel \
    scipy \
    seaborn \
    h5py \
    pyarrow \
    polars \
    pyyaml

# =============================================================================
# liftover Chain Files
# =============================================================================
# Download UCSC liftOver chain files for build conversion
RUN mkdir -p /opt/liftover && \
    wget -q -O /opt/liftover/hg19ToHg38.over.chain.gz \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" && \
    wget -q -O /opt/liftover/hg38ToHg19.over.chain.gz \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# =============================================================================
# EasyColoc Installation
# =============================================================================
WORKDIR /opt

# Copy EasyColoc source
COPY . /opt/EasyColoc
WORKDIR /opt/EasyColoc

# Make CLI entrypoints executable
RUN chmod +x easycoloc tools/*.sh

# =============================================================================
# Default Working Directory
# =============================================================================
WORKDIR /data

# =============================================================================
# Entrypoint
# =============================================================================
ENTRYPOINT ["/opt/EasyColoc/easycoloc"]

# =============================================================================
# Health Check
# =============================================================================
# Verify R and key packages are available
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD R --vanilla -e "library(data.table); library(coloc); library(susieR); cat('OK\n')" || exit 1

# =============================================================================
# Usage Examples
# =============================================================================
#
# Build the image:
#   docker build -t easycoloc:latest .
#
# Run with mounted volumes:
#   docker run -v /host/project:/work -v /host/ref:/ref:ro easycoloc:latest doctor \
#     --global /work/config/global.yml --gwas /work/config/gwas.yml --qtl /work/config/qtl.yml
#
# Run interactively:
#   docker run -it --entrypoint bash -v /host/project:/work easycoloc:latest
#
# Run a managed project:
#   docker run -v /host/project:/work -v /host/ref:/ref:ro easycoloc:latest \
#     run --managed --global /work/config/global.yml --gwas /work/config/gwas.yml --qtl /work/config/qtl.yml --output-dir /work/results
#
# =============================================================================
