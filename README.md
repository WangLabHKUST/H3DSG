# H3DSG: H3 K27-mutant diffuse spinal cord gliomas

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)

## Status
Active Development

## Introduction

This repository contains the code for the multi-omics and spatial single cell investigation of H3-mutant Diffuse Spinal cord Gliomas. We analyzed bulk DNA sequencing, RNA sequencing, DNA methylation array, ATACseq, H3K27me3-CUT&Tag data, and single cell RNA-seq and ATACseq data, as well as CODEX multiplex imaging data. The data were used for multi-omics clustering, survival analysis, cell type annotation, cell of origin analysis and spatial niche identification. We also developed a machine learning classifier based on selected clinical and multi-omics features. Code for these analyses are included in each corresponding folder. A small dataset is also provided to demo the code.

## System requirements and dependencies

Preprocessing of raw sequencing data requires high-performance clusters. We used a cluster equipped with Linux CentOS 7 (kernel version 3.10.0-1062.el7.x86_64) with 40 cores, 256 GB RAM, and at least 100 TB storage.

The remaining code can be run on a desktop. We tested code on a MacBook Pro with macOS Sonoma 14.3, 32 GB RAM and 1 TB storage. We also used QuPath 0.5.0 (https://qupath.github.io/).

Practical environments used for the current repository snapshot:

- R 4.5.1 was used for final runs. Package versions are listed in `env/R_packages.tsv`, and full session info is in `env/R_session.txt`.
- Python 3.12 (Anaconda) package versions are listed in `env/python_packages.txt`.
- `CODEX_Analysis/env/` and `Single_Cell_RNAseq/env/` contain environment notes specific to those analyses.
- `velocity_pipeline_TL6/environment.yml` defines a separate conda environment (Python 3.10, scVelo, velocyto).

## Installation

To install, either download a zip file from https://github.com/WangLabHKUST/H3DSG, or clone via:

```bash
git clone https://github.com/WangLabHKUST/H3DSG.git
```

Install R dependencies according to each analysis folder's scripts (some folders provide `required_packages.R`).
For a repository-level baseline, use:

```r
pkgs <- read.delim("env/R_packages.tsv")$package
install.packages(setdiff(pkgs, rownames(installed.packages())))
# For Bioconductor packages when needed, e.g.:
# BiocManager::install("ComplexHeatmap")
```

Install Python dependencies with:

```bash
pip install -r env/python_requirements.txt
```

## Demo

For each analysis, enter the corresponding directory and run scripts there. Paths inside scripts are generally relative to the script folder.

Example:

```bash
cd WES_Purity_comp
Rscript purity_barplot.R
```

## Raw data availability

The data generated in this study has been deposited in the [Genome Sequence Archive](https://ngdc.cncb.ac.cn/gsa-human/). The accession number for DNA and RNA sequencing data is HRA008075, and the accession number of the DNA methylation data is OMIX006956. The data is not released at the moment but will be released upon publication of the manuscript.

## Contact

For any questions, please contact Professor Jiguang Wang via email: jgwang AT ust DOT hk
