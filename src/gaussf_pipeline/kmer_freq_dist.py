import argparse
import csv
import os
from collections import defaultdict

def process_fasta_in_memory(input_file):
    """
    Reads a FASTA file, extracts gene symbol and transcript ID, formats headers
    as 'gene_symbol|transcript_id' in memory, and returns lists of modified
    headers and corresponding sequences.
    Skips entries missing the 'gene_symbol:' field or the transcript ID.
    """
    modified_headers = []
    sequences = []
    current_sequence = ""
    valid_header_processed = False # Flag to track if the last processed header was valid
    skipped_no_symbol = 0 # Counter for skipped entries due to missing gene symbol

    try:
        with open(input_file, 'r') as infile:
            for line_num, line in enumerate(infile, 1): # Add line number for better error reporting
                line = line.strip()
                if not line: # Skip empty lines
                    continue

                if line.startswith('>'):
                    # Store the previous sequence if it exists and belonged to a valid header
                    if current_sequence and valid_header_processed:
                        sequences.append(current_sequence)

                    current_sequence = "" # Reset sequence accumulator
                    valid_header_processed = False # Reset flag for the new header

                    parts = line.split()
                    try:
                        gene_symbol = None
                        transcript_id = None
                        # Loop through parts to find gene_symbol and transcript_id
                        for part in parts:
                            if part.startswith('gene_symbol:'):
                                # Extract gene symbol
                                gene_symbol_parts = part.split(':', 1)
                                if len(gene_symbol_parts) > 1:
                                    gene_symbol = gene_symbol_parts[1]
                            # The transcript ID is typically the first part after '>'
                            if part.startswith('>'):
                                transcript_id = part[1:] # Store without '>'

                        # --- MODIFICATION START: Check for gene_symbol ---
                        if gene_symbol and transcript_id:
                            # Store the modified header *without* the leading '>'
                            # Format: gene_symbol|transcript_id
                            new_header = f"{gene_symbol}|{transcript_id}"
                            modified_headers.append(new_header)
                            valid_header_processed = True # Mark this header as valid and processed
                        elif transcript_id and not gene_symbol:
                            # Skip if gene_symbol is missing, but log it
                            # print(f"Skipping entry on line {line_num}: Missing 'gene_symbol:' field in header: {line}")
                            skipped_no_symbol += 1
                            valid_header_processed = False # Ensure sequence is skipped
                        else:
                            # Handle other missing info cases (e.g., missing transcript ID)
                            missing_parts = []
                            if not transcript_id: missing_parts.append("transcript ID (first part after '>')")
                            if not gene_symbol: missing_parts.append("'gene_symbol:' field") # Should be caught above, but good practice
                            raise ValueError(f"Header missing required information: {', '.join(missing_parts)}")
                        # --- MODIFICATION END ---

                    except (IndexError, ValueError) as e:
                        print(f"Skipping malformed header on line {line_num}: {line} - Error: {e}")
                        # Keep valid_header_processed as False, sequence lines following this will be ignored

                elif valid_header_processed: # Only add sequence lines if the *last* header processed was valid
                    current_sequence += line

            # Add the last sequence if it exists and belongs to a valid header
            if current_sequence and valid_header_processed:
                sequences.append(current_sequence)

    except FileNotFoundError:
        print(f"Error: Input FASTA file not found at {input_file}")
        return [], [] # Return empty lists on file error
    except Exception as e:
        print(f"An unexpected error occurred during FASTA processing: {e}")
        return [], [] # Return empty lists on other errors

    # Report skipped entries due to missing gene symbol
    if skipped_no_symbol > 0:
        print(f"Note: Skipped {skipped_no_symbol} FASTA entries because they lacked the 'gene_symbol:' field.")


    # Final check for consistency
    if len(modified_headers) != len(sequences):
         # This condition might occur if the file ends abruptly after a header without sequence,
         # or if the last header was invalid/skipped.
         print(f"Warning: Final count mismatch between processed headers ({len(modified_headers)}) and sequences ({len(sequences)}). "
               f"This might be due to the last entry or skipped entries.")
         # Depending on requirements, you might want to adjust lists here, but often proceeding is acceptable.

    return modified_headers, sequences

def sanitize_filename(header):
    """Replaces characters unsuitable for filenames."""
    # Replace pipe, slashes, and backslashes with underscore
    sanitized = header.replace('|', '_').replace('/', '_').replace('\\', '_')
    # Optionally remove or replace other problematic characters like colons, spaces etc.
    # sanitized = sanitized.replace(':', '_').replace(' ', '_')
    return sanitized

# --- Main script logic starts here ---

if __name__ == "__main__":
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a FASTA file for kmer analysis, generating CSV files with kmer counts and transcript occurrences, maintaining kmer order. CSV filenames use gene_symbol_transcriptID.")
    parser.add_argument('--input_fasta', required=True, help="Path to the original input FASTA file.")
    parser.add_argument('--output_dir', required=True, help="Path to the output directory where CSV files will be saved.")
    parser.add_argument('--kmer_length', type=int, default=50, help="Length of the kmers to analyze (default: 50).")
    parser.add_argument('--threshold', type=int, default=3000, help="Maximum number of kmer rows to write to each output CSV file (default: 3000).")
    args = parser.parse_args()

    # --- 1. Process FASTA in memory ---
    print(f"Processing FASTA file: {args.input_fasta}")
    # transcript_headers now contain "gene_symbol|transcript_id"
    transcript_headers, transcripts = process_fasta_in_memory(args.input_fasta)

    if not transcript_headers or not transcripts:
        print("No valid transcript data processed (ensure entries have 'gene_symbol:' field). Exiting.")
        exit(1) # Exit with an error code if FASTA processing failed or yielded no data

    print(f"Successfully processed {len(transcript_headers)} transcripts with gene symbols in memory.")

    # --- 2. Prepare output directory ---
    output_directory = args.output_dir
    try:
        os.makedirs(output_directory, exist_ok=True)
        print(f"Output directory set to: {output_directory}")
    except OSError as e:
        print(f"Error creating output directory {output_directory}: {e}")
        exit(1)

    # --- 3. Kmer analysis ---
    kmer_length = args.kmer_length
    if kmer_length <= 0:
        print("Error: kmer_length must be a positive integer.")
        exit(1)

    global_kmer_counts = defaultdict(int)
    kmer_transcript_sets = defaultdict(set) # Stores indices of transcripts containing the kmer

    print(f"Counting {kmer_length}-mers globally across all processed transcripts...")
    # Kmer counting and tracking transcripts that contain them
    skipped_sequences_count = 0
    for isoform_index, sequence in enumerate(transcripts):
        if len(sequence) < kmer_length:
            # Optionally print a warning for skipped sequences
            # print(f"Warning: Sequence for header '{transcript_headers[isoform_index]}' is shorter than kmer length ({len(sequence)} < {kmer_length}). Skipping kmer analysis for this transcript.")
            skipped_sequences_count += 1
            continue # Skip sequences shorter than kmer length
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            # Basic validation for kmers (optional, e.g., skip kmers with 'N')
            # if 'N' in kmer or 'n' in kmer:
            #     continue
            global_kmer_counts[kmer] += 1
            kmer_transcript_sets[kmer].add(isoform_index) # Store index

    if skipped_sequences_count > 0:
        print(f"Note: Skipped kmer analysis for {skipped_sequences_count} sequences shorter than kmer length ({kmer_length}).")
    print(f"Global kmer counting complete. Found {len(global_kmer_counts)} unique kmers.")

    print(f"Generating CSV files in: {output_directory}")
    # --- Creating CSV files for each isoform ---
    csv_generated_count = 0
    for isoform_index, header in enumerate(transcript_headers):
        # header is now "gene_symbol|transcript_id"
        sequence = transcripts[isoform_index]

        # Skip if the sequence was too short for kmer analysis (already checked, but good practice)
        if len(sequence) < kmer_length:
            continue

        # --- MODIFICATION START: Use sanitized header directly for filename ---
        # sanitize_filename will convert "gene_symbol|transcript_id" to "gene_symbol_transcript_id"
        sanitized_header = sanitize_filename(header)
        output_csv_path = os.path.join(output_directory, sanitized_header + '_kmers.csv')
        # --- MODIFICATION END ---


        # --- Kmer selection and ordering logic for the current isoform ---

        # Pass 1.1: Find the minimum global frequency among kmers present in *this* isoform
        # Pass 1.2: Store kmers from this isoform in their original order along with global freq
        min_global_frequency_for_isoform = float('inf')
        ordered_kmers_in_isoform = [] # List to store (kmer, global_freq) in sequence order
        has_kmers = False
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            # Optional: Skip kmers with 'N' here as well if desired for output
            # if 'N' in kmer or 'n' in kmer:
            #     continue
            global_freq = global_kmer_counts.get(kmer, 0) # Use .get for safety, though kmer should exist
            if global_freq > 0: # Only consider kmers actually counted globally
                 min_global_frequency_for_isoform = min(min_global_frequency_for_isoform, global_freq)
                 ordered_kmers_in_isoform.append((kmer, global_freq))
                 has_kmers = True

        if not has_kmers or min_global_frequency_for_isoform == float('inf'):
            # print(f"No valid {kmer_length}-mers found or processed for {header}. Skipping CSV generation.")
            continue # Skip if no kmers could be generated or none met criteria

        # Pass 2: Collect potential rows (kmers with the minimum global frequency for this isoform),
        # calculate local counts, and format 'Present_in_Transcripts'. Maintain order.
        rows_to_write_ordered = []
        local_kmer_counts_for_isoform = defaultdict(int)

        # Calculate local counts specifically for this isoform
        for kmer, _ in ordered_kmers_in_isoform:
            local_kmer_counts_for_isoform[kmer] += 1

        # Filter based on min_global_frequency and build rows in order
        for kmer, global_freq in ordered_kmers_in_isoform:
            if global_freq == min_global_frequency_for_isoform:
                local_freq = local_kmer_counts_for_isoform[kmer]

                # Format 'Present_in_Transcripts' field using sanitized headers
                transcript_indices = sorted(list(kmer_transcript_sets.get(kmer, set()))) # Get indices, sort them
                # Use the already modified headers (gene_symbol|transcript_id) and sanitize them for the column
                sanitized_transcript_list = [sanitize_filename(transcript_headers[idx]) for idx in transcript_indices]

                if len(sanitized_transcript_list) > 1:
                     transcripts_containing_kmer = '-'.join(sanitized_transcript_list)
                elif sanitized_transcript_list: # Should always be at least one
                     transcripts_containing_kmer = sanitized_transcript_list[0]
                else:
                     # This case should ideally not happen if global_freq > 0
                     transcripts_containing_kmer = "Error_TranscriptNotFound"
                     print(f"Warning: Kmer {kmer} found with global_freq {global_freq} but no transcript index in kmer_transcript_sets for {header}.")

                rows_to_write_ordered.append((kmer, local_freq, global_freq, transcripts_containing_kmer))

        # Pass 3: Filter based on the most frequent 'Present_in_Transcripts' content among the selected kmers
        if rows_to_write_ordered:
            transcript_content_freq = defaultdict(int)
            for _, _, _, transcripts_content in rows_to_write_ordered:
                transcript_content_freq[transcripts_content] += 1

            max_transcript_content_frequency = 0
            if transcript_content_freq: # Check if the dictionary is not empty
                 max_transcript_content_frequency = max(transcript_content_freq.values())

            # Write to CSV only for rows matching the highest "Present_in_Transcripts" frequency,
            # maintaining the original kmer order.
            try:
                with open(output_csv_path, mode='w', newline='') as csv_file:
                    csv_writer = csv.writer(csv_file)
                    csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])

                    written_count = 0
                    for row in rows_to_write_ordered: # Iterate through the ordered list
                        # Check if the 'Present_in_Transcripts' string for this row is one of the most frequent ones
                        if transcript_content_freq[row[3]] == max_transcript_content_frequency:
                            csv_writer.writerow(row)
                            written_count += 1
                            if written_count >= args.threshold:
                                # print(f"Reached threshold ({args.threshold}) for {header}. Stopping CSV writing.")
                                break
                if written_count > 0:
                    csv_generated_count += 1
                # else:
                    # print(f"No kmers met the final filtering criteria for {header}.")

            except IOError as e:
                print(f"Error writing CSV file {output_csv_path}: {e}")
            except Exception as e:
                 print(f"An unexpected error occurred while writing CSV for {header}: {e}")
        # else:
            # print(f"No kmers met the minimum global frequency criteria ({min_global_frequency_for_isoform}) for {header}.")


    print(f"\nKmer analysis complete. Generated {csv_generated_count} CSV files in: {output_directory}")
