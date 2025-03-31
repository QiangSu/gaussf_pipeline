import os
import pytest
from gaussf_pipeline.kmer_freq_dist import main as kmer_freq_dist_main  # Adjust import if needed

def test_kmer_freq_dist(tmp_path):
    input_fasta = "test_data/Homo_sapiens.GRCh38_first_200_transcripts.fa"
    output_dir = tmp_path / "index_output"
    output_dir.mkdir()

    # Run the script
    args = [
        "--input_fasta", str(input_fasta),
        "--output_dir", str(output_dir),
        "--kmer_length", "50",
        "--threshold", "3000"
    ]
    kmer_freq_dist_main(args)

    # Check if output files are created
    assert any(f.endswith(".csv") for f in os.listdir(output_dir)), "No CSV files generated"

