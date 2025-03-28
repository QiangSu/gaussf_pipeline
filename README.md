# GaussF Pipeline

A collection of scripts to perform k-mer counting, normalization, and Gaussian fitting for transcript abundance estimation.

## Description

This pipeline consists of four main steps:
1.  **Index Generation:** Creates k-mer reference files from a transcript FASTA file (`kmer_freq_dist.py`).
2.  **K-mer Counting:** Counts k-mers from FASTQ sequencing data based on the reference index (`kmer_counter.py`).
3.  **Merging & Normalization:** Merges counts back to the reference and calculates normalized counts (`merge_normalize.py`).
4.  **Abundance Estimation:** Applies Gaussian fitting to estimate transcript abundance (`gaussf_tpm.py`).

## Installation

```bash
# Clone the repository (or download)
git clone https://github.com/your_username/gaussf_pipeline.git
cd gaussf_pipeline

# Install using pip (editable mode recommended for development)
pip install -e .

