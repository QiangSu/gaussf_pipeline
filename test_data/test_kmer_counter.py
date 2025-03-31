#!/usr/bin/env python3
import os
from gaussf_pipeline.kmer_counter import main

if __name__ == "__main__":
    # Define test input and output paths
    fastq_path = "test_data/99272-N_extracted_100_raw_2.fastq.gz"
    csv_input_dir = "test_data/csv_input"
    csv_output_dir = "test_data/output/kmer_counter"
    
    # Ensure directories exist
    os.makedirs(csv_input_dir, exist_ok=True)
    os.makedirs(csv_output_dir, exist_ok=True)
    
    # Simulate an input CSV file (required by kmer_counter.py)
    sample_csv = os.path.join(csv_input_dir, "sample_kmers.csv")
    if not os.path.exists(sample_csv):
        with open(sample_csv, "w", newline="") as f:
            import csv
            writer = csv.writer(f)
            writer.writerow(["kmer"])
            writer.writerow(["ATCG" * 12 + "A"])  # Example 50-mer

    # Set command-line arguments matching kmer_counter.py's expectations
    sys.argv = [
        "",  # Script name (placeholder)
        "--fastq_path", fastq_path,
        "--num_threads", "2",  # Minimal threads for testing
        "--chunk_size", "100",  # Small chunk size for testing
        "--csv_input_dir", csv_input_dir,
        "--csv_output_dir", csv_output_dir
    ]
    
    # Run the main function
    main()
