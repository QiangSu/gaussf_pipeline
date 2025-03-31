#!/usr/bin/env python3
import os
from gaussf_pipeline.kmer_freq_dist import main

if __name__ == "__main__":
    input_file = "test_data/Homo_sapiens.GRCh38_first_200_transcripts.fa"
    output_dir = "test_data/output/kmer_freq_dist"
    os.makedirs(output_dir, exist_ok=True)
    # Pass arguments as a list to main()
    main(["--input_fasta", input_file, "--output_dir", output_dir, "--kmer_length", "50", "--threshold", "3000"])
