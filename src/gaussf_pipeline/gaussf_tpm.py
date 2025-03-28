import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import argparse
import re
from tqdm import tqdm  # Import tqdm for the progress bar
import sys # Good practice to import sys for potential exit calls or stderr

# --- Function Definitions ---

# Define a function to judge the suitability of the mean before fitting
def suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
    # Ensure w_fitted_local is positive before comparison
    if w_fitted_local <= 0:
        return False
    return xc_fitted_local > 0.5 * w_fitted_local

# Function to calculate total_normalized_kmer_count across all CSV files
def sum_normalized_kmer_counts(input_directory):
    total_normalized_kmer_count = 0
    print(f"Calculating total normalized k-mer count from: {input_directory}") # Added print for clarity

    try:
        # List all relevant CSV files in the input directory
        csv_files = [f for f in os.listdir(input_directory) if f.endswith('merged_normalized.csv')] # Match the files processed later

        if not csv_files:
            print(f"Warning: No '*_merged_normalized.csv' files found in {input_directory}. Total count will be zero.", file=sys.stderr)
            return 0

        # Loop through each CSV file
        for csv_file in tqdm(csv_files, desc="Summing counts", unit="file", leave=False): # Use tqdm here too
            file_path = os.path.join(input_directory, csv_file)
            try:
                # Read the CSV file into a DataFrame
                df = pd.read_csv(file_path)

                # Check if the required column exists
                if 'Normalized_K-mer_Count' in df.columns:
                    # Sum the 'Normalized_K-mer_Count' column, handling potential NaNs
                    file_normalized_kmer_sum = df['Normalized_K-mer_Count'].sum(skipna=True)
                    # Add the file sum to the total sum
                    total_normalized_kmer_count += file_normalized_kmer_sum
                else:
                    print(f"Warning: 'Normalized_K-mer_Count' column missing in {csv_file}. Skipping its contribution.", file=sys.stderr)

            except pd.errors.EmptyDataError:
                print(f"Warning: File {csv_file} is empty. Skipping.", file=sys.stderr)
            except Exception as e:
                print(f"Warning: Could not process file {csv_file} during sum calculation: {e}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input directory not found: {input_directory}", file=sys.stderr)
        sys.exit(1) # Exit if the input directory doesn't exist
    except Exception as e:
        print(f"An unexpected error occurred during sum calculation: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Total Normalized K-mer Count calculated: {total_normalized_kmer_count}")
    # Handle division by zero case later if needed
    if total_normalized_kmer_count == 0:
         print("Warning: Total Normalized K-mer Count is zero. TPM-like normalization might result in division by zero or NaN.", file=sys.stderr)

    return total_normalized_kmer_count

# GC content calculation function
def calculate_gc_content(kmer):
    # Handle potential non-string input or empty strings
    if not isinstance(kmer, str) or not kmer:
        return 0.0
    gc_count = kmer.upper().count('G') + kmer.upper().count('C') # Use upper for robustness
    total_bases = len(kmer)
    if total_bases == 0:
        return 0.0
    gc_content_percent = (gc_count / total_bases) * 100
    return round(gc_content_percent, 2)

# Gaussian CDF definition
def gaussian_cdf(x, A0, A, xc, w):
    # Add safeguard against non-positive w
    if w <= 0:
        # Return a large value to discourage the fitter from choosing w <= 0
        return 1e10 * (1 + abs(x)) # Or np.inf if numpy is imported
    return A0 + A * norm.cdf(x, loc=xc, scale=w) # Use norm.cdf directly

# Gaussian CDF with fixed parameters for Normalized K-mer Count
def gaussian_cdf_fixed(x, A0, A, xc_fixed, w_fixed):
    # w_fixed should already be validated before being passed here
    return A0 + A * norm.cdf(x, loc=xc_fixed, scale=w_fixed)

# Gaussian CDF with fixed parameters for Count
def gaussian_cdf_fixed_count(x, A0, A, xc_fixed, w_fixed):
    # w_fixed should already be validated before being passed here
    return A0 + A * norm.cdf(x, loc=xc_fixed, scale=w_fixed)

# Extract gene name and transcript ID from filename
def extract_gene_transcript_id(filename):
    # Match for mouse (Mus) data with transcript ID and hyphen in gene name
    match_mus = re.search(r'([\w.-]+)_(ENSMUST[0-9]+(?:.\d+)?)_kmers', filename) # Allow dots/hyphens in gene, capture version with transcript
    # Match for human (Homo sapiens) data with transcript ID and hyphen in gene name
    match_human = re.search(r'([\w.-]+)_(ENST[0-9]+(?:.\d+)?)_kmers', filename) # Allow dots/hyphens in gene, capture version with transcript
    # Match for files with no transcript ID (only gene name with hyphen/dot)
    match_no_transcript = re.search(r'([\w.-]+)_kmers', filename)

    if match_mus:
        gene_name = match_mus.group(1)
        transcript_id = match_mus.group(2)
        return gene_name, transcript_id
    elif match_human:
        gene_name = match_human.group(1)
        transcript_id = match_human.group(2)
        return gene_name, transcript_id
    elif match_no_transcript:
        # Check if it didn't match the transcript patterns first
        if not (re.search(r'_(ENSMUST[0-9]+)', filename) or re.search(r'_(ENST[0-9]+)', filename)):
             gene_name = match_no_transcript.group(1)
             return gene_name, "-" # Return hyphen for missing transcript ID
    # Fallback if no pattern matches
    base = os.path.basename(filename)
    name_part = base.split('_kmers')[0] # Try a simple split
    print(f"Warning: Could not parse standard gene/transcript from '{filename}'. Using '{name_part}' as Gene_Name.", file=sys.stderr)
    return name_part, "Unknown"


# --- Main Execution Logic ---
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Analyze GC content, fit Gaussian CDF, and estimate abundance.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the '*_merged_normalized.csv' files.")
    parser.add_argument('--output', type=str, required=True, help="Path and name of the output CSV file to save the results.")
    parser.add_argument('--threshold', type=int, default=10, help="Minimum number of distinct GC content data points required for fitting. Default is 10.")
    args = parser.parse_args()

    # --- Start of Execution ---
    print("Starting Gaussian fitting analysis...")

    # Calculate total_normalized_kmer_count across all relevant CSV files in the input directory
    # This needs to be done *before* the loop if normalization factor is global
    total_normalized_kmer_count = sum_normalized_kmer_counts(args.input)

    # List to store the results
    results = []

    # Get the list of files to process
    try:
        all_files_in_dir = os.listdir(args.input)
    except FileNotFoundError:
        print(f"Error: Input directory not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    files_to_process = sorted([f for f in all_files_in_dir if f.endswith("merged_normalized.csv")])

    if not files_to_process:
        print(f"Error: No '*_merged_normalized.csv' files found in the input directory: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files_to_process)} files to process.")

    # Loop through each file in the directory with a progress bar
    for filename in tqdm(files_to_process, desc="Processing files", unit="file"):
        filepath = os.path.join(args.input, filename)
        gene_name, transcript_id = "Unknown", "Unknown" # Initialize defaults

        try:
            # Extract gene name and transcript ID from the filename
            gene_name, transcript_id = extract_gene_transcript_id(filename)

            # Read the CSV file into a DataFrame
            df = pd.read_csv(filepath)

            # Basic checks on the DataFrame
            if df.empty:
                print(f"Warning: File {filename} is empty. Skipping.", file=sys.stderr)
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                    'Global_Frequency': 'N/A', 'Present_in_Transcripts': 'N/A', 'Transcript_Length': 'N/A',
                    'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                    'Sum or Fitted A (Abundance) for Count': '0.00',
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': 'Empty input file'
                })
                continue

            required_columns = ['kmer', 'Transcript_Length', 'Local_Frequency', 'Normalized_K-mer_Count', 'Count', 'Global_Frequency', 'Present_in_Transcripts']
            if not all(col in df.columns for col in required_columns):
                missing_cols = [col for col in required_columns if col not in df.columns]
                print(f"Warning: File {filename} is missing required columns: {missing_cols}. Skipping.", file=sys.stderr)
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                    'Global_Frequency': df.get('Global_Frequency', pd.Series(['N/A']))[0], # Try to get some info
                    'Present_in_Transcripts': df.get('Present_in_Transcripts', pd.Series(['N/A']))[0],
                    'Transcript_Length': df.get('Transcript_Length', pd.Series(['N/A']))[0],
                    'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                    'Sum or Fitted A (Abundance) for Count': '0.00',
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': f'Missing columns: {missing_cols}'
                })
                continue

            # Extract metadata (assuming it's constant per file)
            transcript_length = df['Transcript_Length'].iloc[0]
            global_frequency = df['Global_Frequency'].iloc[0]
            present_in_transcripts = df['Present_in_Transcripts'].iloc[0]

            # Calculate and add GC content to the DataFrame
            df['GC_Content'] = df['kmer'].apply(calculate_gc_content)

            # Group by GC content and sum frequencies
            gc_content_data = df.groupby('GC_Content').agg({
                'Local_Frequency': 'sum',
                'Normalized_K-mer_Count': 'sum',
                'Count': 'sum'
            }).reset_index()

            # Filter out rows where GC_Content might be NaN or Inf if calculation failed
            gc_content_data = gc_content_data.dropna(subset=['GC_Content'])
            gc_content_data = gc_content_data[gc_content_data['GC_Content'].apply(lambda x: isinstance(x, (int, float)))]


            # Check if there are at least args.threshold distinct GC contents *after* grouping
            if len(gc_content_data) < args.threshold:
                report_reason = 'Not enough distinct GC contents'
                sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
                # Normalize the sum by multiplying by 1000000/total_normalized_kmer_count (handle division by zero)
                normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count else 0
                sum_count = gc_content_data['Count'].sum()
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                    'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': report_reason
                })
                continue # Skip to the next file

            # --- Proceed with Fitting ---
            # Calculate cumulative sums
            gc_content_data_sorted = gc_content_data.sort_values(by='GC_Content')
            # Ensure cumulative sums are calculated correctly even if counts are zero
            gc_content_data_sorted['Cumulative_Local_Frequency'] = gc_content_data_sorted['Local_Frequency'].cumsum()
            gc_content_data_sorted['Cumulative_Normalized_Count'] = gc_content_data_sorted['Normalized_K-mer_Count'].cumsum()
            gc_content_data_sorted['Cumulative_Count'] = gc_content_data_sorted['Count'].cumsum()

            # Get the data for fitting
            x_data = gc_content_data_sorted['GC_Content'].values # Use .values for numpy arrays
            y_data_local = gc_content_data_sorted['Cumulative_Local_Frequency'].values
            y_data_normalized = gc_content_data_sorted['Cumulative_Normalized_Count'].values
            y_data_count = gc_content_data_sorted['Cumulative_Count'].values

            # Check for sufficient variance in x_data
            if x_data.std() < 1e-6:
                 report_reason = 'GC content variance too low for fitting'
                 sum_normalized_kmer_count = gc_content_data_sorted['Normalized_K-mer_Count'].sum()
                 normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count else 0
                 sum_count = gc_content_data_sorted['Count'].sum()
                 results.append({
                     'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                     'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                     'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                     'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                     'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                     'Report': report_reason
                 })
                 continue

            # Fit the Gaussian CDF to Cumulative Local Frequency to get initial xc and w
            # Provide reasonable bounds, especially for w > 0
            bounds_local = ([min(y_data_local)-abs(min(y_data_local)), 0, min(x_data), 1e-3], # A0_min, A_min=0, xc_min, w_min > 0
                            [max(y_data_local)+abs(max(y_data_local)), (max(y_data_local) - min(y_data_local)) * 1.5, max(x_data), max(x_data)-min(x_data)]) # A0_max, A_max, xc_max, w_max
            initial_guesses_local = [min(y_data_local), max(y_data_local) - min(y_data_local), x_data.mean(), x_data.std() if x_data.std() > 1e-3 else 1.0]
            # Ensure initial w guess is positive
            if initial_guesses_local[3] <= 0: initial_guesses_local[3] = 1.0

            try:
                popt_local, pcov_local = curve_fit(gaussian_cdf, x_data, y_data_local, p0=initial_guesses_local, bounds=bounds_local, maxfev=5000)
                A0_fitted_local, A_fitted_local, xc_fitted_local, w_fitted_local = popt_local

                # Check if the fitted mean and standard deviation meet the suitability criteria
                if suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
                    # Additional fitting for Normalized K-mer Count and Count using fixed xc and w
                    # Bounds for A0 and A
                    bounds_norm = ([min(y_data_normalized)-abs(min(y_data_normalized)), 0],
                                   [max(y_data_normalized)+abs(max(y_data_normalized)), (max(y_data_normalized) - min(y_data_normalized)) * 1.5])
                    initial_guesses_normalized = [min(y_data_normalized), max(y_data_normalized) - min(y_data_normalized)]
                    popt_normalized, pcov_normalized = curve_fit(
                        lambda x, A0, A: gaussian_cdf_fixed(x, A0, A, xc_fitted_local, w_fitted_local),
                        x_data, y_data_normalized, p0=initial_guesses_normalized, bounds=bounds_norm, maxfev=5000
                    )
                    A0_fitted_normalized, A_fitted_normalized = popt_normalized

                    # Normalize A_fitted_normalized (handle division by zero)
                    A_fitted_normalized_tpm = (A_fitted_normalized * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count else 0

                    bounds_count = ([min(y_data_count)-abs(min(y_data_count)), 0],
                                    [max(y_data_count)+abs(max(y_data_count)), (max(y_data_count) - min(y_data_count)) * 1.5])
                    initial_guesses_count = [min(y_data_count), max(y_data_count) - min(y_data_count)]
                    popt_count, pcov_count = curve_fit(
                        lambda x, A0, A: gaussian_cdf_fixed_count(x, A0, A, xc_fitted_local, w_fitted_local),
                        x_data, y_data_count, p0=initial_guesses_count, bounds=bounds_count, maxfev=5000
                    )
                    A0_fitted_count, A_fitted_count = popt_count

                    # Append successful fitting results
                    results.append({
                        'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                        'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                        'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(A_fitted_normalized_tpm),
                        'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(A_fitted_count),
                        'Fixed Mean (xc)': '{:.2f}'.format(xc_fitted_local),
                        'Fixed Standard Deviation (w)': '{:.2f}'.format(w_fitted_local),
                        'Report': 'OK'
                    })
                else:
                    # Use the specific reason from the criteria function if possible, otherwise generic message
                    report_reason = f'Unsuitable Fit Parameters (xc={xc_fitted_local:.2f}, w={w_fitted_local:.2f})'
                    raise ValueError(report_reason) # Raise error to be caught below

            except (RuntimeError, ValueError) as e:
                # Handle fitting failures or unsuitable parameters by reporting sums
                error_message = str(e)
                sum_normalized_kmer_count = gc_content_data_sorted['Normalized_K-mer_Count'].sum()
                normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count else 0
                sum_count = gc_content_data_sorted['Count'].sum()
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                    'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': f'Fit Failed or Unsuitable - {error_message}'
                })

        except FileNotFoundError:
            print(f"Error: File not found {filepath}. Skipping.", file=sys.stderr)
            # Append minimal info if possible
            results.append({
                'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                'Global_Frequency': 'N/A', 'Present_in_Transcripts': 'N/A', 'Transcript_Length': 'N/A',
                'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                'Sum or Fitted A (Abundance) for Count': '0.00',
                'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                'Report': 'File not found during processing'
            })
        except pd.errors.EmptyDataError:
             print(f"Warning: File {filename} is empty. Skipping.", file=sys.stderr)
             results.append({
                 'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                 'Global_Frequency': 'N/A', 'Present_in_Transcripts': 'N/A', 'Transcript_Length': 'N/A',
                 'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                 'Sum or Fitted A (Abundance) for Count': '0.00',
                 'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                 'Report': 'Empty input file'
             })
        except Exception as e:
            print(f"Error processing file {filename}: {e}", file=sys.stderr)
            # Append error information
            results.append({
                'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id,
                'Global_Frequency': 'Error', 'Present_in_Transcripts': 'Error', 'Transcript_Length': 'Error',
                'Sum or Fitted A (Abundance) for Normalized Count': 'Error',
                'Sum or Fitted A (Abundance) for Count': 'Error',
                'Fixed Mean (xc)': 'Error', 'Fixed Standard Deviation (w)': 'Error',
                'Report': f'Unexpected Error: {e}'
            })

    # --- Final Output ---
    if not results:
        print("No results were generated. Check input files and logs.", file=sys.stderr)
        sys.exit(1) # Exit if no results at all

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Ensure the output directory exists
    try:
        output_directory = os.path.dirname(args.output)
        if output_directory: # Only create if path includes a directory
            os.makedirs(output_directory, exist_ok=True)
        # Save to CSV file specified by the command-line argument
        results_df.to_csv(args.output, index=False)
        print(f"\nResults successfully saved to {args.output}")
    except OSError as e:
        print(f"Error: Could not create output directory or save file at {args.output}: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to save results to {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

# --- Call the main function only when the script is executed directly ---
if __name__ == "__main__":
    main()

