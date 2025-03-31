import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import argparse
import re
from tqdm import tqdm  # Import tqdm for the progress bar
import sys # Good practice to import sys for potential exit calls or stderr
# import numpy as np # Import if using np.inf in gaussian_cdf

# --- Function Definitions ---

# Define a function to judge the suitability of the mean before fitting
def suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
    # Ensure w_fitted_local is positive before comparison
    if w_fitted_local <= 0:
        return False
    # Basic sanity check: mean should be reasonably larger than half the width,
    # especially relevant if GC content is expected to be > 0.
    # This helps avoid fits where the peak is near or below zero GC content.
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
        # Using np.inf requires importing numpy as np
        # return np.inf
        # Ensure x is treated as float for arithmetic operation
        try:
            x_float = float(x)
        except (ValueError, TypeError):
             # Handle cases where x might not be convertible to float, though unlikely if data is cleaned
             return 1e10 # Return a large constant if conversion fails
        return 1e10 * (1 + abs(x_float)) # Avoid inf to prevent potential issues downstream
    return A0 + A * norm.cdf(x, loc=xc, scale=w) # Use norm.cdf directly

# Gaussian CDF with fixed parameters for Normalized K-mer Count
def gaussian_cdf_fixed(x, A0, A, xc_fixed, w_fixed):
    # w_fixed should already be validated before being passed here
    # No need for w <= 0 check here as it's checked before fixing
    return A0 + A * norm.cdf(x, loc=xc_fixed, scale=w_fixed)

# Gaussian CDF with fixed parameters for Count
def gaussian_cdf_fixed_count(x, A0, A, xc_fixed, w_fixed):
    # w_fixed should already be validated before being passed here
    # No need for w <= 0 check here as it's checked before fixing
    return A0 + A * norm.cdf(x, loc=xc_fixed, scale=w_fixed)

# --- NEW Function to Extract Gene Name and Transcript ID by underscore position ---
def extract_gene_transcript_id(filename):
    """
    Extracts Gene Name and Transcript ID based on underscore positions.
    Gene Name: Content before the first underscore.
    Transcript ID: Content strictly between the first and second underscore.

    Args:
        filename (str): The input filename (just the name, not the path).

    Returns:
        tuple: (gene_name, transcript_id)
               Returns (FileNamePart, "-") if fewer than two underscores exist.
    """
    # Split the filename by the underscore character
    parts = filename.split('_')
    num_parts = len(parts)

    if num_parts == 1:
        # No underscores found. Treat the whole relevant part as the gene name.
        # Attempt to remove common suffixes for a cleaner name.
        gene_name = filename
        if gene_name.endswith('_merged_normalized.csv'):
             gene_name = gene_name[:-len('_merged_normalized.csv')]
        elif gene_name.endswith('_kmers.csv'):
             gene_name = gene_name[:-len('_kmers.csv')]
        elif gene_name.endswith('.csv'):
             gene_name = gene_name[:-len('.csv')]
        transcript_id = "-" # No underscores means no transcript ID by this rule
    elif num_parts == 2:
        # Exactly one underscore. Gene name is the first part.
        # No content *between* first and second underscore.
        gene_name = parts[0]
        transcript_id = "-"
    else: # num_parts >= 3
        # Two or more underscores. Gene name is the first part.
        # Transcript ID is the second part (between first and second underscore).
        gene_name = parts[0]
        transcript_id = parts[1]

    return gene_name, transcript_id
# --- End of NEW Function ---

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
            # --- Use the NEW extraction function ---
            gene_name, transcript_id = extract_gene_transcript_id(filename)
            # ---

            # Read the CSV file into a DataFrame
            df = pd.read_csv(filepath)

            # Basic checks on the DataFrame
            if df.empty:
                print(f"Warning: File {filename} is empty. Skipping.", file=sys.stderr)
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
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
                # Attempt to get metadata even if some columns are missing
                global_freq_val = df.get('Global_Frequency', pd.Series(['N/A'])).iloc[0] if not df.empty else 'N/A'
                present_in_val = df.get('Present_in_Transcripts', pd.Series(['N/A'])).iloc[0] if not df.empty else 'N/A'
                length_val = df.get('Transcript_Length', pd.Series(['N/A'])).iloc[0] if not df.empty else 'N/A'

                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
                    'Global_Frequency': global_freq_val,
                    'Present_in_Transcripts': present_in_val,
                    'Transcript_Length': length_val,
                    'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                    'Sum or Fitted A (Abundance) for Count': '0.00',
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': f'Missing columns: {missing_cols}'
                })
                continue

            # Extract metadata (assuming it's constant per file)
            # Use .iloc[0] safely after checking df is not empty and columns exist
            transcript_length = df['Transcript_Length'].iloc[0]
            global_frequency = df['Global_Frequency'].iloc[0]
            present_in_transcripts = df['Present_in_Transcripts'].iloc[0]

            # Calculate and add GC content to the DataFrame
            df['GC_Content'] = df['kmer'].apply(calculate_gc_content)

            # Group by GC content and sum frequencies
            gc_content_data = df.groupby('GC_Content', dropna=False).agg({ # Include potential NaN group initially
                'Local_Frequency': 'sum',
                'Normalized_K-mer_Count': 'sum',
                'Count': 'sum'
            }).reset_index()

            # Filter out rows where GC_Content might be NaN or Inf if calculation failed or kmer was invalid
            gc_content_data = gc_content_data.dropna(subset=['GC_Content'])
            # Ensure GC_Content is numeric after potential grouping issues
            gc_content_data = gc_content_data[pd.to_numeric(gc_content_data['GC_Content'], errors='coerce').notna()]
            gc_content_data['GC_Content'] = gc_content_data['GC_Content'].astype(float)


            # Check if there are at least args.threshold distinct GC contents *after* grouping and cleaning
            if len(gc_content_data) < args.threshold:
                report_reason = f'Not enough distinct GC contents ({len(gc_content_data)} < {args.threshold})'
                sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
                # Normalize the sum by multiplying by 1000000/total_normalized_kmer_count (handle division by zero)
                normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count > 0 else 0.0
                sum_count = gc_content_data['Count'].sum()
                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
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

            # Check for sufficient variance in x_data and non-empty data
            if len(x_data) < 2 or x_data.std() < 1e-6:
                 report_reason = 'GC content variance too low or insufficient points for fitting'
                 sum_normalized_kmer_count = gc_content_data_sorted['Normalized_K-mer_Count'].sum()
                 normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count > 0 else 0.0
                 sum_count = gc_content_data_sorted['Count'].sum()
                 results.append({
                     'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
                     'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                     'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                     'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                     'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                     'Report': report_reason
                 })
                 continue

            # Fit the Gaussian CDF to Cumulative Local Frequency to get initial xc and w
            # Provide reasonable bounds, especially for w > 0 and A > 0
            # Calculate dynamic bounds based on data range
            y_local_min, y_local_max = y_data_local.min(), y_data_local.max()
            y_local_range = y_local_max - y_local_min if y_local_max > y_local_min else 1.0 # Avoid zero range
            x_min, x_max = x_data.min(), x_data.max()
            x_range = x_max - x_min if x_max > x_min else 1.0 # Avoid zero range

            bounds_local = ([y_local_min - y_local_range * 0.5, 0, x_min, 1e-3], # A0_min, A_min=0, xc_min, w_min > 0
                            [y_local_max + y_local_range * 0.5, y_local_range * 1.5, x_max, x_range]) # A0_max, A_max, xc_max, w_max

            # Sensible initial guesses
            initial_guesses_local = [y_local_min, y_local_range, x_data.mean(), x_data.std() if x_data.std() > 1e-3 else 1.0]
            # Ensure initial w guess is positive and within bounds
            initial_guesses_local[3] = max(bounds_local[0][3], min(bounds_local[1][3], initial_guesses_local[3]))
            # Ensure initial xc guess is within bounds
            initial_guesses_local[2] = max(bounds_local[0][2], min(bounds_local[1][2], initial_guesses_local[2]))


            try:
                popt_local, pcov_local = curve_fit(gaussian_cdf, x_data, y_data_local, p0=initial_guesses_local, bounds=bounds_local, maxfev=10000, method='trf') # Try 'trf' or 'dogbox' for bounds
                A0_fitted_local, A_fitted_local, xc_fitted_local, w_fitted_local = popt_local

                # Check if the fitted mean and standard deviation meet the suitability criteria
                if suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
                    # Additional fitting for Normalized K-mer Count and Count using fixed xc and w

                    # --- Normalized Count Fit ---
                    y_norm_min, y_norm_max = y_data_normalized.min(), y_data_normalized.max()
                    y_norm_range = y_norm_max - y_norm_min if y_norm_max > y_norm_min else 1e-9 # Avoid zero, allow small values

                    bounds_norm = ([y_norm_min - y_norm_range * 0.5, 0], # A0_min, A_min=0
                                   [y_norm_max + y_norm_range * 0.5, y_norm_range * 1.5]) # A0_max, A_max
                    initial_guesses_normalized = [y_norm_min, y_norm_range]
                    initial_guesses_normalized[1] = max(bounds_norm[0][1], min(bounds_norm[1][1], initial_guesses_normalized[1])) # Ensure A guess > 0

                    popt_normalized, pcov_normalized = curve_fit(
                        lambda x, A0, A: gaussian_cdf_fixed(x, A0, A, xc_fixed=xc_fitted_local, w_fixed=w_fitted_local),
                        x_data, y_data_normalized, p0=initial_guesses_normalized, bounds=bounds_norm, maxfev=10000, method='trf'
                    )
                    A0_fitted_normalized, A_fitted_normalized = popt_normalized

                    # Normalize A_fitted_normalized (handle division by zero)
                    A_fitted_normalized_tpm = (A_fitted_normalized * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count > 0 else 0.0

                    # --- Count Fit ---
                    y_count_min, y_count_max = y_data_count.min(), y_data_count.max()
                    y_count_range = y_count_max - y_count_min if y_count_max > y_count_min else 1.0 # Avoid zero

                    bounds_count = ([y_count_min - y_count_range * 0.5, 0], # A0_min, A_min=0
                                    [y_count_max + y_count_range * 0.5, y_count_range * 1.5]) # A0_max, A_max
                    initial_guesses_count = [y_count_min, y_count_range]
                    initial_guesses_count[1] = max(bounds_count[0][1], min(bounds_count[1][1], initial_guesses_count[1])) # Ensure A guess > 0

                    popt_count, pcov_count = curve_fit(
                        lambda x, A0, A: gaussian_cdf_fixed_count(x, A0, A, xc_fixed=xc_fitted_local, w_fixed=w_fitted_local),
                        x_data, y_data_count, p0=initial_guesses_count, bounds=bounds_count, maxfev=10000, method='trf'
                    )
                    A0_fitted_count, A_fitted_count = popt_count

                    # Append successful fitting results
                    results.append({
                        'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
                        'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                        'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(A_fitted_normalized_tpm),
                        'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(A_fitted_count),
                        'Fixed Mean (xc)': '{:.2f}'.format(xc_fitted_local),
                        'Fixed Standard Deviation (w)': '{:.2f}'.format(w_fitted_local),
                        'Report': 'OK'
                    })
                else:
                    # Fit deemed unsuitable based on xc vs w relationship
                    report_reason = f'Unsuitable Fit Parameters (xc={xc_fitted_local:.2f}, w={w_fitted_local:.2f}) fails xc > 0.5*w check'
                    # Raise ValueError to be caught by the except block below, forcing sum calculation
                    raise ValueError(report_reason)

            except (RuntimeError, ValueError) as e:
                # Handle fitting failures (RuntimeError from curve_fit) or unsuitable parameters (ValueError raised above)
                # Report sums as fallback
                error_message = str(e)
                # Check if we have the sorted data, otherwise calculate sum from original grouped data
                if 'gc_content_data_sorted' in locals():
                    sum_normalized_kmer_count = gc_content_data_sorted['Normalized_K-mer_Count'].sum()
                    sum_count = gc_content_data_sorted['Count'].sum()
                else:
                    sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
                    sum_count = gc_content_data['Count'].sum()

                normalized_sum = (sum_normalized_kmer_count * 1_000_000 / total_normalized_kmer_count) if total_normalized_kmer_count > 0 else 0.0

                # Determine the report message based on the error
                report_prefix = 'Fit Failed' if isinstance(e, RuntimeError) else 'Fit Unsuitable'
                final_report = f'{report_prefix} - Using Sum - {error_message}'

                results.append({
                    'File': filename, 'Gene_Name': gene_name, 'Transcript_ID': transcript_id, # Use parsed names
                    'Global_Frequency': global_frequency, 'Present_in_Transcripts': present_in_transcripts, 'Transcript_Length': transcript_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                    'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                    'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                    'Report': final_report
                })

        except FileNotFoundError:
            print(f"Error: File not found {filepath}. Skipping.", file=sys.stderr)
            # Use names parsed by the new function even if file not found later
            gene_name_parsed, transcript_id_parsed = extract_gene_transcript_id(filename)
            results.append({
                'File': filename, 'Gene_Name': gene_name_parsed, 'Transcript_ID': transcript_id_parsed,
                'Global_Frequency': 'N/A', 'Present_in_Transcripts': 'N/A', 'Transcript_Length': 'N/A',
                'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                'Sum or Fitted A (Abundance) for Count': '0.00',
                'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                'Report': 'File not found during processing loop'
            })
        except pd.errors.EmptyDataError:
             print(f"Warning: File {filename} is empty. Skipping.", file=sys.stderr)
             # Use names parsed by the new function
             gene_name_parsed, transcript_id_parsed = extract_gene_transcript_id(filename)
             results.append({
                 'File': filename, 'Gene_Name': gene_name_parsed, 'Transcript_ID': transcript_id_parsed,
                 'Global_Frequency': 'N/A', 'Present_in_Transcripts': 'N/A', 'Transcript_Length': 'N/A',
                 'Sum or Fitted A (Abundance) for Normalized Count': '0.00',
                 'Sum or Fitted A (Abundance) for Count': '0.00',
                 'Fixed Mean (xc)': 'N/A', 'Fixed Standard Deviation (w)': 'N/A',
                 'Report': 'Empty input file (pd.errors.EmptyDataError)'
             })
        except Exception as e:
            print(f"Error processing file {filename}: {type(e).__name__} - {e}", file=sys.stderr)
            # Use names parsed by the new function
            gene_name_parsed, transcript_id_parsed = extract_gene_transcript_id(filename)
            # Attempt to retrieve metadata if possible, otherwise use 'Error'
            try:
                 global_freq_val = df.get('Global_Frequency', pd.Series(['Error'])).iloc[0] if 'df' in locals() and not df.empty else 'Error'
                 present_in_val = df.get('Present_in_Transcripts', pd.Series(['Error'])).iloc[0] if 'df' in locals() and not df.empty else 'Error'
                 length_val = df.get('Transcript_Length', pd.Series(['Error'])).iloc[0] if 'df' in locals() and not df.empty else 'Error'
            except:
                 global_freq_val, present_in_val, length_val = 'Error', 'Error', 'Error'

            results.append({
                'File': filename, 'Gene_Name': gene_name_parsed, 'Transcript_ID': transcript_id_parsed,
                'Global_Frequency': global_freq_val,
                'Present_in_Transcripts': present_in_val,
                'Transcript_Length': length_val,
                'Sum or Fitted A (Abundance) for Normalized Count': 'Error',
                'Sum or Fitted A (Abundance) for Count': 'Error',
                'Fixed Mean (xc)': 'Error', 'Fixed Standard Deviation (w)': 'Error',
                'Report': f'Unexpected Error: {type(e).__name__} - {e}'
            })

    # --- Final Output ---
    if not results:
        print("No results were generated. Check input files and logs.", file=sys.stderr)
        sys.exit(1) # Exit if no results at all

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Define column order for clarity
    column_order = [
        'File', 'Gene_Name', 'Transcript_ID', 'Transcript_Length',
        'Global_Frequency', 'Present_in_Transcripts',
        'Sum or Fitted A (Abundance) for Normalized Count',
        'Sum or Fitted A (Abundance) for Count',
        'Fixed Mean (xc)', 'Fixed Standard Deviation (w)',
        'Report'
    ]
    # Reorder columns, handling potential missing columns if errors occurred
    results_df = results_df.reindex(columns=column_order, fill_value='N/A')


    # Ensure the output directory exists
    try:
        output_directory = os.path.dirname(args.output)
        if output_directory and not os.path.exists(output_directory): # Only create if path includes a directory and it doesn't exist
             print(f"Creating output directory: {output_directory}")
             os.makedirs(output_directory, exist_ok=True)
        elif not output_directory:
             print("Output path does not specify a directory. Saving to current directory.")

        # Save to CSV file specified by the command-line argument
        results_df.to_csv(args.output, index=False)
        print(f"\nResults successfully saved to {args.output}")
    except OSError as e:
        print(f"Error: Could not create output directory or save file at {args.output}: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to save results to {args.output}: {type(e).__name__} - {e}", file=sys.stderr)
        sys.exit(1)

# --- Call the main function only when the script is executed directly ---
if __name__ == "__main__":
    main()
