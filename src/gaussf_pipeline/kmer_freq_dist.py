import argparse
import csv
import os
from collections import defaultdict

def process_fasta_in_memory(input_file):
    """
    Reads a FASTA file, extracts gene symbol and transcript ID, formats headers
    as 'gene_symbol|transcript_id' in memory, and returns lists of modified
    headers and sequences.
    """
    modified_headers = []
    sequences = []
    current_sequence = []
    current_header = None

    try:
        with open(input_file, 'r') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_sequence:
                        sequences.append(''.join(current_sequence))
                        current_sequence = []
                    header_parts = line[1:].split()
                    transcript_id = header_parts[0]  # First part is transcript ID
                    gene_symbol = None
                    for part in header_parts:
                        if part.startswith('gene_symbol:'):
                            gene_symbol = part.split(':')[1]
                            break
                    if gene_symbol:
                        current_header = f"{gene_symbol}|{transcript_id}"
                        modified_headers.append(current_header)
                    else:
                        current_header = None  # Skip entries without gene_symbol
                elif current_header:
                    current_sequence.append(line)
            if current_header and current_sequence:
                sequences.append(''.join(current_sequence))
    except IOError as e:
        print(f"Error reading FASTA file {input_file}: {e}")
        exit(1)

    return modified_headers, sequences

def sanitize_filename(header):
    """Replace problematic characters in header for safe filenames."""
    return header.replace('|', '_').replace(':', '_').replace('/', '_').replace('\\', '_')

def main(args=None):
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Process a FASTA file for kmer analysis, generating CSV files with kmer counts and transcript occurrences, maintaining kmer order. CSV filenames use gene_symbol_transcriptID.")
    parser.add_argument('--input_fasta', required=True, help="Path to the original input FASTA file.")
    parser.add_argument('--output_dir', required=True, help="Path to the output directory where CSV files will be saved.")
    parser.add_argument('--kmer_length', type=int, default=50, help="Length of the kmers to analyze (default: 50).")
    parser.add_argument('--threshold', type=int, default=3000, help="Maximum number of kmer rows to write to each output CSV file (default: 3000).")
    args = parser.parse_args(args if args is not None else None)  # Allow passing args as a list

    # --- 1. Process FASTA in memory ---
    print(f"Processing FASTA file: {args.input_fasta}")
    transcript_headers, transcripts = process_fasta_in_memory(args.input_fasta)

    if not transcript_headers or not transcripts:
        print("No valid transcript data processed (ensure entries have 'gene_symbol:' field). Exiting.")
        exit(1)

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
    kmer_transcript_sets = defaultdict(set)

    print(f"Counting {kmer_length}-mers globally across all processed transcripts...")
    skipped_sequences_count = 0
    for isoform_index, sequence in enumerate(transcripts):
        if len(sequence) < kmer_length:
            skipped_sequences_count += 1
            continue
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            global_kmer_counts[kmer] += 1
            kmer_transcript_sets[kmer].add(isoform_index)

    if skipped_sequences_count > 0:
        print(f"Note: Skipped kmer analysis for {skipped_sequences_count} sequences shorter than kmer length ({kmer_length}).")
    print(f"Global kmer counting complete. Found {len(global_kmer_counts)} unique kmers.")

    print(f"Generating CSV files in: {output_directory}")
    csv_generated_count = 0
    for isoform_index, header in enumerate(transcript_headers):
        sequence = transcripts[isoform_index]
        if len(sequence) < kmer_length:
            continue

        sanitized_header = sanitize_filename(header)
        output_csv_path = os.path.join(output_directory, sanitized_header + '_kmers.csv')

        min_global_frequency_for_isoform = float('inf')
        ordered_kmers_in_isoform = []
        has_kmers = False
        for i in range(len(sequence) - kmer_length + 1):
            kmer = sequence[i:i + kmer_length]
            global_freq = global_kmer_counts.get(kmer, 0)
            if global_freq > 0:
                min_global_frequency_for_isoform = min(min_global_frequency_for_isoform, global_freq)
                ordered_kmers_in_isoform.append((kmer, global_freq))
                has_kmers = True

        if not has_kmers or min_global_frequency_for_isoform == float('inf'):
            continue

        rows_to_write_ordered = []
        local_kmer_counts_for_isoform = defaultdict(int)
        for kmer, _ in ordered_kmers_in_isoform:
            local_kmer_counts_for_isoform[kmer] += 1

        for kmer, global_freq in ordered_kmers_in_isoform:
            if global_freq == min_global_frequency_for_isoform:
                local_freq = local_kmer_counts_for_isoform[kmer]
                transcript_indices = sorted(list(kmer_transcript_sets.get(kmer, set())))
                sanitized_transcript_list = [sanitize_filename(transcript_headers[idx]) for idx in transcript_indices]
                transcripts_containing_kmer = '-'.join(sanitized_transcript_list) if len(sanitized_transcript_list) > 1 else sanitized_transcript_list[0] if sanitized_transcript_list else "Error_TranscriptNotFound"
                if transcripts_containing_kmer == "Error_TranscriptNotFound":
                    print(f"Warning: Kmer {kmer} found with global_freq {global_freq} but no transcript index in kmer_transcript_sets for {header}.")
                rows_to_write_ordered.append((kmer, local_freq, global_freq, transcripts_containing_kmer))

        if rows_to_write_ordered:
            transcript_content_freq = defaultdict(int)
            for _, _, _, transcripts_content in rows_to_write_ordered:
                transcript_content_freq[transcripts_content] += 1

            max_transcript_content_frequency = max(transcript_content_freq.values()) if transcript_content_freq else 0

            try:
                with open(output_csv_path, mode='w', newline='') as csv_file:
                    csv_writer = csv.writer(csv_file)
                    csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])
                    written_count = 0
                    for row in rows_to_write_ordered:
                        if transcript_content_freq[row[3]] == max_transcript_content_frequency:
                            csv_writer.writerow(row)
                            written_count += 1
                            if written_count >= args.threshold:
                                break
                if written_count > 0:
                    csv_generated_count += 1
            except IOError as e:
                print(f"Error writing CSV file {output_csv_path}: {e}")

    print(f"\nKmer analysis complete. Generated {csv_generated_count} CSV files in: {output_directory}")

if __name__ == "__main__":
    main()
