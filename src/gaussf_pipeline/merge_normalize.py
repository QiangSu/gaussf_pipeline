import pandas as pd
import os
import glob
import argparse # Keep argparse import
import sys # Import sys for potential error exits
from tqdm import tqdm

# Updated normalization function (removed total_read_count)
def normalize_frequency(kmer_count, total_csv_kmers, read_length, k):
    # Basic check to prevent division by zero or negative values if inputs are unexpected
    denominator = total_csv_kmers * (read_length - k + 1)
    if denominator <= 0:
        # Decide how to handle this: return 0, raise error, or print warning
        print(f"Warning: Invalid denominator components: total_csv_kmers={total_csv_kmers}, read_length={read_length}, k={k}", file=sys.stderr)
        return 0
    multiplier = 2 * 1000 / denominator
    return kmer_count * multiplier

# Function to merge k-mer count files to original k-mer CSV and add normalized count
def merge_kmer_counts_to_kmers(kmers_filepath, kmer_counts_filepath, total_csv_kmers, read_length, k, output_filepath):
    try:
        # Load the original kmer data and kmer counts
        kmers_df = pd.read_csv(kmers_filepath)
        # Check if kmer counts file exists before trying to read
        if not os.path.exists(kmer_counts_filepath):
            print(f"Warning: K-mer counts file not found, skipping: {kmer_counts_filepath}", file=sys.stderr)
            return # Skip processing this file pair
        kmer_counts_df = pd.read_csv(kmer_counts_filepath)

        # Check for empty dataframes
        if kmers_df.empty or kmer_counts_df.empty:
            print(f"Warning: One or both input CSVs are empty, skipping merge for: {os.path.basename(kmers_filepath)}", file=sys.stderr)
            return

        # Ensure required columns exist
        if 'kmer' not in kmers_df.columns:
             print(f"Error: 'kmer' column missing in {kmers_filepath}", file=sys.stderr)
             return
        if 'K-mer' not in kmer_counts_df.columns or 'Count' not in kmer_counts_df.columns:
             print(f"Error: 'K-mer' or 'Count' column missing in {kmer_counts_filepath}", file=sys.stderr)
             return

        # Merge the dataframes on the 'kmer' column
        # Use how='left' to keep all kmers from the reference, even if not found in counts
        merged_df = kmers_df.merge(kmer_counts_df, left_on='kmer', right_on='K-mer', how='left')

        # Drop the redundant 'K-mer' column from the merged dataframe if it exists
        if 'K-mer' in merged_df.columns:
            merged_df.drop('K-mer', axis=1, inplace=True)

        # Fill NaN counts with 0 if a kmer from reference wasn't found in counts
        merged_df['Count'].fillna(0, inplace=True)
        # Ensure Count is integer type after fillna
        merged_df['Count'] = merged_df['Count'].astype(int)

        # Add normalized k-mer counts
        merged_df['Normalized_K-mer_Count'] = merged_df['Count'].apply(normalize_frequency, args=(total_csv_kmers, read_length, k))

        # Add the 'Transcript_Length' column (total kmers in the reference file)
        merged_df['Transcript_Length'] = total_csv_kmers

        # Write the merged data to a new CSV file
        merged_df.to_csv(output_filepath, index=False)

    except pd.errors.EmptyDataError:
        print(f"Warning: Empty CSV file encountered, skipping: {kmers_filepath} or {kmer_counts_filepath}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: File not found during processing: {kmers_filepath} or {kmer_counts_filepath}", file=sys.stderr)
    except Exception as e:
        print(f"Error processing {os.path.basename(kmers_filepath)}: {e}", file=sys.stderr)


# Main function - NOW CONTAINS ARGUMENT PARSING
def main():
    # --- Argument parsing moved inside main ---
    parser = argparse.ArgumentParser(description='Merge k-mer counts into original k-mer reference CSV files and add normalized counts.')
    parser.add_argument('--kmer_reference_directory', type=str, required=True, help='Directory containing the k-mer reference CSV files (*_kmers.csv).')
    parser.add_argument('--kmer_counts_directory', type=str, required=True, help='Directory containing the k-mer counts CSV files (*_kmers_counts.csv).')
    parser.add_argument('--output_directory', type=str, required=True, help='Output directory for storing merged CSV files (*_merged_normalized.csv).')
    parser.add_argument('--read_length', type=int, default=150, help='Read length of FASTQ sequences.')
    parser.add_argument('--k', type=int, default=50, help='K-mer length.')
    args = parser.parse_args() # Parse args from sys.argv

    # Use parsed arguments directly via args object
    kmer_reference_directory = args.kmer_reference_directory
    kmer_counts_directory = args.kmer_counts_directory
    output_directory = args.output_directory
    read_length = args.read_length
    k = args.k
    # --- End of argument parsing section ---

    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True) # Use exist_ok=True

    # Process all the files
    kmer_files_pattern = os.path.join(kmer_reference_directory, "*_kmers.csv")
    kmer_files = glob.glob(kmer_files_pattern)

    if not kmer_files:
        print(f"Warning: No '*_kmers.csv' files found in {kmer_reference_directory}", file=sys.stderr)
        return # Exit if no files to process

    print(f"Found {len(kmer_files)} k-mer reference files in {kmer_reference_directory}")

    # Use tqdm for progress bar
    for kmer_filepath in tqdm(kmer_files, desc="Merging and Normalizing"):
        try:
            # Calculate the number of k-mers in the current reference CSV file
            # Add check for empty file before calculating shape
            if os.path.getsize(kmer_filepath) == 0:
                 print(f"Warning: Skipping empty k-mer reference file: {kmer_filepath}", file=sys.stderr)
                 continue
            total_csv_kmers = pd.read_csv(kmer_filepath).shape[0]
            if total_csv_kmers == 0:
                 print(f"Warning: K-mer reference file has 0 rows: {kmer_filepath}", file=sys.stderr)
                 # Decide if you want to skip or proceed (proceeding might lead to division by zero in normalize_frequency)
                 continue # Skip if no kmers

            # Compute the filename of the corresponding kmer counts file
            # Assumes counts file name is derived from reference file name
            basename = os.path.basename(kmer_filepath) # e.g., "transcript_xyz_kmers.csv"
            ref_name_part = basename.replace('_kmers.csv', '') # e.g., "transcript_xyz"
            # --- Adjust this line based on the EXACT naming convention of your count files ---
            # Example: If count files are named "transcript_xyz_kmers_counts.csv"
            kmer_counts_filename = f"{ref_name_part}_kmers_counts.csv"
            # Example: If count files are named "transcript_xyz_counts.csv"
            # kmer_counts_filename = f"{ref_name_part}_counts.csv"
            # --- End adjustment ---
            kmer_counts_filepath = os.path.join(kmer_counts_directory, kmer_counts_filename)

            # Define the output file path in the output directory
            output_filename = f"{ref_name_part}_merged_normalized.csv"
            output_filepath = os.path.join(output_directory, output_filename)

            # Merge the kmer counts into the kmer file and add normalized counts
            merge_kmer_counts_to_kmers(kmer_filepath, kmer_counts_filepath, total_csv_kmers, read_length, k, output_filepath)

        except pd.errors.EmptyDataError:
             print(f"Warning: Skipping empty k-mer reference file after initial check: {kmer_filepath}", file=sys.stderr)
        except Exception as e:
             print(f"Error processing file {kmer_filepath}: {e}", file=sys.stderr)

    print(f"\nProcessing complete. Merged files saved to: {output_directory}")


# Run script - NOW SIMPLY CALLS main()
if __name__ == "__main__":
    main() # Call main without arguments

