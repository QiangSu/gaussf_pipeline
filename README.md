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
git clone https://github.com/QiangSu/gaussf_pipeline.git
cd gaussf_pipeline

# Install using pip (editable mode recommended for development)
pip install -e .

# 1. Generate reference index
# (Command will likely change to something like 'gaussf_index_gen --input ...')
python src/gaussf_pipeline/kmer_freq_dist.py --input_fasta /path/to/transcriptome.fa --output_dir /path/to/index_output --kmer_length 50 --threshold 3000

# 2. Count k-mers from FASTQ
# (Command will likely change to something like 'gaussf_kmer_count --fastq_path ...')
python src/gaussf_pipeline/kmer_counter.py --fastq_path /path/to/reads.fastq.gz --num_threads 30 --chunk_size 100000 --csv_input_dir /path/to/index_output --csv_output_dir /path/to/counts_output

# 3. Merge and Normalize
# (Command will likely change to something like 'gaussf_merge --kmer_reference_directory ...')
python src/gaussf_pipeline/merge_normalize.py --kmer_reference_directory /path/to/index_output --kmer_counts_directory /path/to/counts_output --output_directory /path/to/merged_output --read_length 150 --k 50

# 4. Estimate Abundance
# (Command will likely change to something like 'gaussf_abundance --input ...')
python src/gaussf_pipeline/gaussf_tpm.py --threshold 10 --input /path/to/merged_output --output /path/to/results.csv
