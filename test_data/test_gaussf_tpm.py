import os
import pytest
from gaussf_pipeline.gaussf_tpm import main as gaussf_tpm_main  # Adjust import if needed

def test_gaussf_tpm(tmp_path):
    input_dir = tmp_path / "merged_output"
    input_dir.mkdir()
    output_file = tmp_path / "results.csv"

    # Mock input file
    with open(input_dir / "mock_merged.csv", "w") as f:
        f.write("transcript,count,length\nT1,100,500\n")

    # Run the script
    args = [
        "--threshold", "10",
        "--input", str(input_dir),
        "--output", str(output_file)
    ]
    gaussf_tpm_main(args)

    # Check if output file exists and has content
    assert output_file.exists(), "Results file not created"
    assert output_file.stat().st_size > 0, "Results file is empty"

