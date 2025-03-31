import os
import pytest
from gaussf_pipeline.merge_normalize import main as merge_normalize_main  # Adjust import if needed

def test_merge_normalize(tmp_path):
    kmer_ref_dir = tmp_path / "index_output"
    kmer_ref_dir.mkdir()
    kmer_counts_dir = tmp_path / "counts_output"
    kmer_counts_dir.mkdir()
    output_dir = tmp_path / "merged_output"
    output_dir.mkdir()

    # Mock input files
    with open(kmer_ref_dir / "ref.csv", "w") as f:
        f.write("kmer,count\nATCG,10\n")
    with open(kmer_counts_dir / "counts.csv", "w") as f:
        f.write("kmer,count\nATCG,5\n")

    # Run the script
    args = [
        "--kmer_reference_directory", str(kmer_ref_dir),
        "--kmer_counts_directory", str(kmer_counts_dir),
        "--output_directory", str(output_dir),
        "--read_length", "150",
        "--k", "50"
    ]
    merge_normalize_main(args)

    # Check if output files are created
    assert len(os.listdir(output_dir)) > 0, "No merged output files generated"

