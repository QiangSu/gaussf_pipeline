#!/usr/bin/env python3
import os
import sys
import csv
from gaussf_pipeline.merge_normalize import main

if __name__ == "__main__":
    # Define test directories
    kmer_reference_dir = "test_data/output/kmer_freq_dist"
    kmer_counts_dir = "test_data/output/kmer_counter"
    output_dir = "test_data/output/merge_normalize"
    
    # Ensure directories exist
    os.makedirs(kmer_reference_dir, exist_ok=True)
    os.makedirs(kmer_counts_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Create sample kmer reference file (from kmer_freq_dist.py)
    sample_kmer_file = os.path.join(kmer_reference_dir, "test_transcript_kmers.csv")
    if not os.path.exists(sample_kmer_file):
        with open(sample_kmer_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["kmer"])
            writer.writerow(["ATCG" * 12 + "A"])  # 50-mer example
    
    # Create sample kmer counts file (from kmer_counter.py)
    sample_counts_file = os.path.join(kmer_counts_dir, "test_transcript_kmers_counts.csv")
    if not os.path.exists(sample_counts_file):
        with open(sample_counts_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["K-mer", "Count"])
            writer.writerow(["ATCG" * 12 + "A", 5])  # Example count
    
    # Set command-line arguments matching merge_normalize.py
    sys.argv = [
        "",  # Script name
        "--kmer_reference_directory", kmer_reference_dir,
        "--kmer_counts_directory", kmer_counts_dir,
        "--output_directory", output_dir,
        "--read_length", "150",  # Default value
        "--k", "50"  # Default value
    ]
    
    # Run the main function
    main()
