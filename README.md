# H3DSG: H3 K27-mutant diffuse spinal cord gliomas

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)

## Status
Active Development

## Introduction

This repository contains the code for the multi-omics and spatial single cell investigation of H3-mutant Diffuse Spinal cord Gliomas. We analyzed bulk DNA sequencing, RNA sequencing, DNA methylation array, ATACseq, H3K27me3-CUT&Tag data, and single cell RNA-seq and ATACseq data, as well as CODEX multiplex imaging data. The data were used for multi-omics clustering, survival analysis, cell type annotation, cell of origin analysis and spatial niche identification. We also developed a machine learning classifier based on selected clinical and multi-omics features. Code for these analysis were included in each corresponding folder. A small dataset is also provided to demo the code.

## System requirements and dependencies
Preprocessing of the raw sequencing data requires high performance clusters. We used a cluster equipped with Linux CentOS 7 (kernel version 3.10.0-1062.el7.x86_64) with 40 cores, 256 GB RAM and at least 100TB storage. Code for these preprocessing steps are in the "Preprocessing" folder.

The other parts of the code can be run on a desktop. We tested the code on a MacBook Pro with macOS Sonoma 14.3, 32 GB RAM and 1 TB storage. Required dependencies include R 4.4.1 and python 3.9.6. The required R and python packages are included in each of the code snippets. We also used QuPath 0.5.0 which is downloadable at https://qupath.github.io/. 

## Installation

To install, the user can either download a zip file from https://github.com/WangLabHKUST/H3DSG, or by cloning the repository via
```
git clone https://github.com/WangLabHKUST/H3DSG.git
```
To install the dependent R packages, please refer to the manual of each individual package.

## Demo
For each analysis, the user can enter the directory and run the R code inside the folder. The code will read the demo data provided in the folder, run the analysis, and generate result files within the folder.

## Raw data availability
The data generated in this study has been deposited in the [Genome Sequence Archive](https://ngdc.cncb.ac.cn/gsa-human/). The accession number for DNA and RNA sequencing data is HRA008075, and the accession number of the DNA methylation data is OMIX006956.

## Contact
For any questions, please contact Professor Jiguang Wang via email: jgwang AT ust DOT hk
