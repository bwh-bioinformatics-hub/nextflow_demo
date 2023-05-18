# nextflow_demo
Demo of a Nextflow Pipeline for bulk-RNA Seq

**************************

# 1. Introduction

### Overview:
Note: This pipeline is designed to be a demo to learn how to construct as well as run a nextflow pipeline.

This pipeline manages a RNA-Seq workflow starting from raw fastq files and converting
them to standard file formats for use by downstream tools. The steps involved are:

* FastQC: Import of data from BAM, SAM or FastQ files (any variant). This tool provides a quick overview to tell you in which areas there may be problems. Creates a summary graphs and tables to quickly assess your data.
* Trim Galore: Adapter Trimming
* STAR: Alignment of reads to reference genome.

<a id="dependencies"></a>

## Dependencies    
This repository uses Nextflow for pipeline managment, Conda for Environment Management, FastQC, STAR, and TrimGalore.
```
We have handled all software requirements using Conda.
The only thing user will need to install is miniconda, which can be installed following this manual:
https://educe-ubc.github.io/conda.html
```
Clone this repo:
```
git clone https://github.com/bwh-bioinformatics-hub/nextflow_demo.git
```
Once you have cloned the repo, to create conda environment with dependencies installed
```
conda env create -f env/environment.yml 
```
Once the Conda environment has been installed to run:
```
nextflow run rna_seq_pipeline.nf with-conda ~/miniconda3/envs/rna_seq_pipeline

