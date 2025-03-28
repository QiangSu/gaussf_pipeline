import gzip
import concurrent.futures
from collections import Counter
from queue import Queue
import csv
import os
import argparse
import gc
from tqdm import tqdm


def producer(file_path, queue, chunk_size):
    # Determine if the file is compressed based on the file extension
    if file_path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(file_path, 'rt') as f:
        sequences = []
        line_count = sum(1 for line in f)  # Total lines in the file for progress bar
        f.seek(0)  # Reset file pointer to the beginning
        pbar = tqdm(total=line_count // 4, desc='Reading sequences')
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            _ = f.readline().strip()  # plus line (ignore)
            _ = f.readline().strip()  # quality line (ignore)
            sequences.append(sequence)
            pbar.update(1)
            if len(sequences) == chunk_size:
                queue.put(sequences)
                sequences = []
                gc.collect()  # Trigger garbage collection
        if sequences:
            queue.put(sequences)
        queue.put(None)  # Signal the end of file to consumers
        pbar.close()


# Consumer function for processing sequences into kmers and counting them
def consumer(queue, k=50):
    local_freq_dict = Counter()
    while True:
        sequences = queue.get()
        if sequences is None:
            queue.put(None)  # Signal next consumer the end of processing
            break
        for sequence in sequences:
            if len(sequence) >= k:  # Safety check for short sequences
                kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
                local_freq_dict.update(kmers)
        del sequences  # Free memory immediately after processing
        gc.collect()  # Trigger garbage collection
    return local_freq_dict


# Function to combine all steps and orchestrate the multiprocessing
def process_and_count(file_path, chunk_size, num_workers):
    queue = Queue(maxsize=num_workers)  # Limit the queue size to avoid excessive memory usage
    global_freq_dict = Counter()

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers + 1) as executor:
        # Setup all consumers
        futures = [executor.submit(consumer, queue, 50) for _ in range(num_workers)]
        # Setup producer
        producer_future = executor.submit(producer, file_path, queue, chunk_size)

        # Wait for the producer to finish (not necessary but keeps things clean)
        concurrent.futures.wait([producer_future])

        # Collect results from consumers
        for fut in concurrent.futures.as_completed(futures):
            global_freq_dict.update(fut.result())

    return global_freq_dict


# Function to read k-mers from a CSV file and write their counts to a new CSV file
def match_kmers_with_counts(global_freq_dict, input_csv_file, output_csv_file):
    with open(input_csv_file, mode='r') as infile:
        reader = csv.DictReader(infile)
        kmers_list = [row['kmer'] for row in reader]

    with open(output_csv_file, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['K-mer', 'Count'])
        for kmer in kmers_list:
            count = global_freq_dict.get(kmer, 0)
            writer.writerow([kmer, count])


# Main function to set up argument parsing
def main():
    parser = argparse.ArgumentParser(description='Process FASTQ file and match k-mers from CSV files.')
    parser.add_argument('--fastq_path', type=str, required=True, help='Path to the FASTQ.gz file')
    parser.add_argument('--num_threads', type=int, required=True, help='Number of worker threads')
    parser.add_argument('--chunk_size', type=int, required=True, help='Chunk size for processing')
    parser.add_argument('--csv_input_dir', type=str, required=True, help='Directory containing input CSV files')
    parser.add_argument('--csv_output_dir', type=str, required=True, help='Directory to save output CSV files')

    args = parser.parse_args()

    # Process and count k-mers
    global_freq_dict = process_and_count(args.fastq_path, args.chunk_size, args.num_threads)

    # Ensure the output directory exists
    os.makedirs(args.csv_output_dir, exist_ok=True)

    # Process each CSV file in the input directory
    input_files = [f for f in os.listdir(args.csv_input_dir) if f.endswith('_kmers.csv')]
    pbar = tqdm(total=len(input_files), desc='Processing CSV files')

    for input_csv_file in input_files:
        input_csv_path = os.path.join(args.csv_input_dir, input_csv_file)
        output_csv_file = input_csv_file.replace('.csv', '_counts.csv')
        output_csv_path = os.path.join(args.csv_output_dir, output_csv_file)
        match_kmers_with_counts(global_freq_dict, input_csv_path, output_csv_path)
        pbar.update(1)

    pbar.close()


if __name__ == '__main__':
    main()
