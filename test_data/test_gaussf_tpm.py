#!/usr/bin/env python3
import os
import sys
import csv
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
from gaussf_pipeline.gaussf_tpm import main

if __name__ == "__main__":
    # Define directories and files
    input_dir = "test_data/output/merge_normalize"
    output_file = "test_data/output/results.csv"
    threshold = "10"

    # Ensure directories exist
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Create a sample input file that matches merge_normalize.py output
    sample_input_file = os.path.join(input_dir, "test_transcript_merged_normalized.csv")
    if not os.path.exists(sample_input_file):
        # Sample data mimicking merge_normalize.py output
        data = {
            "kmer": ["ATCG" * 12 + "A"],
            "Transcript_Length": [150],
            "Local_Frequency": [5.0],
            "Normalized_K-mer_Count": [99.01],  # From previous normalization example
            "Count": [5],
            "Global_Frequency": [0.001],
            "Present_in_Transcripts": [1]
        }
        df = pd.DataFrame(data)
        df.to_csv(sample_input_file, index=False)
        print(f"Created sample input file: {sample_input_file}")

    # Set sys.argv to match gaussf_tpm.py expectations
    sys.argv = [
        "", "--threshold", threshold,
        "--input", input_dir,
        "--output", output_file
    ]

    # Run the main function
    main()
