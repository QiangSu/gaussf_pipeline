import os
import pytest
from gaussf_pipeline.kmer_counter import main as kmer_counter_main  # Adjust import if needed

def test_kmer_counter(tmp_path):
    fastq_path = "test_data/99272-N_extracted_100_raw_2.fastq.gz"
    csv_input_dir = tmp_path / "index_output"
    csv_input_dir.mkdir()
    csv_output_dir = tmp_path / "counts_output"
    csv_output_dir.mkdir()

    # Mock a CSV file in csv_input_dir (simplified for testing)
    with open(csv_input_dir / "mock.csv", "w") as f:
        f.write("kmer,count\nATCG,10\n")

    # Run the script
    args = [
        "--fastq_path", str(fastq_path),
        "--num_threads", "2",  # Reduced for CI
        "--chunk_size", "1000",
        "--csv_input_dir", str(csv_input_dir),
        "--csv_output_dir", str(csv_output_dir)
    ]
    kmer_counter_main(args)

    # Check if output files are created
    assert len(os.listdir(csv_output_dir)) > 0, "No output files generated"

