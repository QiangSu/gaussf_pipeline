# --- 0. Load Necessary Libraries ---
# Ensure you have polyester and Biostrings installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("polyester")
# BiocManager::install("Biostrings")

library(polyester)
library(Biostrings)
# library(readr) # Using base R's read.csv instead

# --- Register Temporary File Cleanup ---
# To ensure temporary files are cleaned up even if the script errors out
temp_files_to_clean <- character(0)
on.exit({
    if (length(temp_files_to_clean) > 0) {
        cleanup_paths <- temp_files_to_clean[file.exists(temp_files_to_clean)]
        if(length(cleanup_paths) > 0) {
             message("Cleaning up temporary files: ", paste(cleanup_paths, collapse=", "))
             unlink(cleanup_paths)
        }
    }
}, add = TRUE)


# --- 1. Define Input and Output Variables ---
reference_transcriptome <- "/home/data/reference_isoform/Homo_sapiens.GRCh38.cdna.all.fa"
# This file only provides the LIST of transcripts to get high base abundance (20).
abundance_list_file <- "/home/data/reference_isoform/99272N_results_file_GC.csv"
output_dir <- "simulated_data_isoform_ASSIGNEDCOUNTS_10000000" # Updated dir name

# --- 2. Define Simulation Parameters ---
# These might need adjustment based on the original experiment or desired simulation setup
set.seed(123) # Set seed for reproducibility
READ_LENGTH <- 150 # Adjust: common read lengths are 75, 100, 125, 150
PAIRED_END <- TRUE # Set to FALSE if you want single-end reads
FRAGMENT_MEAN_LENGTH <- 350 # Adjust: typical mean fragment length (insert size + 2*read_length)
FRAGMENT_SD_LENGTH <- 25 # Adjust: typical standard deviation of fragment length
NUM_REPLICATES_PER_GROUP <- 2 # How many simulated samples (replicates) per group?
NUM_GROUPS <- 1 # Number of experimental groups (e.g., control, treatment). Set > 1 for DE simulation.
TOTAL_SAMPLES <- NUM_GROUPS * NUM_REPLICATES_PER_GROUP
TARGET_TOTAL_READS_PER_SAMPLE <- 10000000 # Desired *approximate* total reads per simulated sample

# *** Define Assigned Base Abundance Levels (BEFORE SCALING) ***
ABUNDANCE_LISTED_TX <- 20 # Base abundance for transcripts listed in abundance_list_file
ABUNDANCE_UNLISTED_TX <- 5  # Base abundance for transcripts from reference NOT listed in abundance_list_file

# --- 3. Read Transcriptome Sequences (ALL OF THEM) ---
message("Reading reference transcriptome FASTA (all transcripts)...")
fasta_sequences <- tryCatch({
    readDNAStringSet(reference_transcriptome)
}, error = function(e) {
    stop("Failed to read reference transcriptome FASTA file '", reference_transcriptome, "': ", e$message)
})
message("Read ", length(fasta_sequences), " sequences from reference FASTA.")
# Extract transcript IDs from FASTA headers (often the first word)
fasta_ids <- names(fasta_sequences)
fasta_ids_simple <- sapply(strsplit(fasta_ids, "\\s+"), `[`, 1) # Get ID before any space
names(fasta_sequences) <- fasta_ids_simple # Update names to simplified IDs
# Remove any sequences with NA or empty names after simplification
valid_fasta_indices <- !is.na(names(fasta_sequences)) & names(fasta_sequences) != ""
fasta_sequences <- fasta_sequences[valid_fasta_indices]
fasta_ids_simple_valid <- names(fasta_sequences) # These are IDs of ALL valid transcripts from reference
message("Using simplified transcript IDs. Kept ", length(fasta_sequences), " sequences with valid IDs. Example: ", head(fasta_ids_simple_valid, 1))

# --- 4. Read List of Transcripts for High Abundance ---
message("Reading list of high-abundance transcripts from: ", abundance_list_file)
message("NOTE: Only the transcript IDs in the first column will be used. Abundance values in this file are ignored.")
abundance_list_data <- tryCatch({
    read.csv(abundance_list_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    stop("Failed to read transcript list CSV file '", abundance_list_file, "': ", e$message)
})

if (ncol(abundance_list_data) < 1) {
    stop("Transcript list file '", abundance_list_file, "' does not have at least one column.")
}
transcript_id_col <- colnames(abundance_list_data)[1]
message("Using column '", transcript_id_col, "' for transcript IDs to assign high base abundance (", ABUNDANCE_LISTED_TX, ").")

listed_tx_ids <- abundance_list_data[[transcript_id_col]]

# Ensure IDs are character and remove NAs/duplicates/empty
if (!is.character(listed_tx_ids)) {
    listed_tx_ids <- as.character(listed_tx_ids)
}
listed_tx_ids <- unique(listed_tx_ids[!is.na(listed_tx_ids) & listed_tx_ids != ""]) # Keep unique, non-NA, non-empty IDs
message("Read ", length(listed_tx_ids), " unique, non-empty transcript IDs from the list file.")

# --- 5. Assign Base Abundances (to ALL transcripts from FASTA) ---
message("Assigning base abundances to ALL ", length(fasta_ids_simple_valid), " transcripts from reference FASTA...")

# Identify which FASTA transcripts are in the high-abundance list
ids_in_list_and_fasta <- intersect(fasta_ids_simple_valid, listed_tx_ids)
num_high_abundance <- length(ids_in_list_and_fasta)

# Identify FASTA transcripts NOT in the list (these get low abundance)
ids_not_in_list <- setdiff(fasta_ids_simple_valid, listed_tx_ids) # These are the "unmatched" transcripts
num_low_abundance <- length(ids_not_in_list)

message("Found ", num_high_abundance, " transcripts from the list file that are also in the reference FASTA (will get base abundance ", ABUNDANCE_LISTED_TX, ").")
message(num_low_abundance, " transcripts from the reference FASTA were NOT in the list file (will get base abundance ", ABUNDANCE_UNLISTED_TX, ").")
num_list_not_in_fasta <- length(setdiff(listed_tx_ids, fasta_ids_simple_valid))
if (num_list_not_in_fasta > 0) {
    message("Note: ", num_list_not_in_fasta, " IDs from the list file were not found in the reference FASTA and will be ignored for abundance assignment.")
}

if (length(fasta_ids_simple_valid) == 0) {
    stop("No valid transcript IDs found in the reference FASTA after initial processing.")
}

# Create a named vector to hold the base counts for ALL valid FASTA transcripts
# Initialize ALL with the low abundance value
base_counts_vector <- setNames(rep(ABUNDANCE_UNLISTED_TX, length(fasta_ids_simple_valid)),
                               fasta_ids_simple_valid)

# Update the counts for the transcripts that were in the list
if (num_high_abundance > 0) {
    base_counts_vector[ids_in_list_and_fasta] <- ABUNDANCE_LISTED_TX
}

message("Assigned base abundance ", ABUNDANCE_LISTED_TX, " to ", sum(base_counts_vector == ABUNDANCE_LISTED_TX), " transcripts.")
message("Assigned base abundance ", ABUNDANCE_UNLISTED_TX, " to ", sum(base_counts_vector == ABUNDANCE_UNLISTED_TX), " transcripts.")

# --- 6. Prepare Count Matrix for Polyester (Using ALL transcripts) ---
message("Preparing count matrix for simulation (includes ALL reference transcripts)...")

# Calculate the total counts based on the assigned base abundances across ALL transcripts
current_total_counts <- sum(base_counts_vector)
message("Sum of assigned base counts across all ", length(base_counts_vector), " transcripts: ", round(current_total_counts))

# Scale these base counts to achieve the target total reads per sample
if (current_total_counts > 0 && TARGET_TOTAL_READS_PER_SAMPLE > 0) {
  scaling_factor <- TARGET_TOTAL_READS_PER_SAMPLE / current_total_counts
  message("Scaling ALL base counts by factor: ", round(scaling_factor, 4), " to target ~", TARGET_TOTAL_READS_PER_SAMPLE, " reads per sample.")
  scaled_counts_vector <- round(base_counts_vector * scaling_factor)
} else {
  message("Warning: Sum of base counts is zero or target reads is zero. Using base counts directly for the matrix.")
  scaled_counts_vector <- base_counts_vector
}

# Ensure counts are non-negative integers after scaling and rounding
scaled_counts_vector[scaled_counts_vector < 0] <- 0
scaled_counts_vector <- round(scaled_counts_vector) # Ensure integer counts

# Create the matrix: rows are ALL transcripts, columns are samples
count_matrix <- matrix(scaled_counts_vector,
                       nrow = length(scaled_counts_vector),
                       ncol = TOTAL_SAMPLES)

rownames(count_matrix) <- names(scaled_counts_vector) # Should be fasta_ids_simple_valid
colnames(count_matrix) <- paste0("sample_", 1:TOTAL_SAMPLES)

message("Created count matrix with dimensions: ", nrow(count_matrix), " transcripts x ", ncol(count_matrix), " samples.")
if(any(is.na(count_matrix))) {
   stop("NA values detected in the final count matrix before simulation. Check base count assignment and scaling.")
}
message("Approximate expected total reads per sample after scaling (colSums): ", paste(round(colSums(count_matrix)), collapse=", "))

zero_count_transcripts <- sum(rowSums(count_matrix) == 0)
if (zero_count_transcripts > 0) {
    message("Note: ", zero_count_transcripts, " transcripts have zero expected counts across all samples after scaling (due to low base abundance and scaling/rounding).")
}
if (nrow(count_matrix) > 0 && nrow(count_matrix) == zero_count_transcripts){
    warning("All transcripts have zero expected counts after scaling. Simulation will likely produce empty files or error.")
}

# --- 7. Prepare FASTA Subset (containing ALL transcripts included in the count matrix) ---
message("Selecting FASTA sequences corresponding to the count matrix (ALL ", nrow(count_matrix), " transcripts)...")
# We select the sequences from the original (filtered) fasta_sequences using the rownames of the final matrix.
fasta_subset <- fasta_sequences[rownames(count_matrix)]

# Double check alignment
if(length(fasta_subset) != nrow(count_matrix)) {
    stop("Fatal Error: Mismatch between number of sequences in final FASTA subset (", length(fasta_subset),
         ") and number of rows in count matrix (", nrow(count_matrix), ").")
}
if(!identical(names(fasta_subset), rownames(count_matrix))) {
    stop("Fatal Error: Names/order mismatch between final FASTA subset and count matrix rownames.")
}
message("Prepared FASTA subset containing ", length(fasta_subset), " sequences for simulation.")


# --- 8. Run Polyester Simulation (Using ALL transcripts) ---

message("Checking inputs before simulation...")
if (!inherits(fasta_subset, "DNAStringSet") || length(fasta_subset) == 0) {
     stop("fasta_subset is invalid or empty before simulation.")
}
if (any(is.na(names(fasta_subset)))) {
    stop("NA names detected in fasta_subset before simulation!")
}
if (nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
    stop("Count matrix is empty. Cannot run simulation. Check input files and matching steps.")
}

message("Writing FASTA subset for simulation (all ", length(fasta_subset), " transcripts) to a temporary file...")
temp_fasta_path <- tempfile(pattern="polyester_fasta_subset_", fileext=".fa")
temp_files_to_clean <- c(temp_files_to_clean, temp_fasta_path) # Register for cleanup
tryCatch({
    writeXStringSet(fasta_subset, filepath=temp_fasta_path, format="fasta")
    message("Temporary FASTA written successfully to: ", temp_fasta_path)
}, error = function(e) {
    stop("Failed to write temporary FASTA file: ", e$message)
})

# Check if the temp file exists and is not empty
if (!file.exists(temp_fasta_path) || file.info(temp_fasta_path)$size == 0) {
    size_info <- if(file.exists(temp_fasta_path)) paste0("size ", file.info(temp_fasta_path)$size, " bytes") else "does not exist"
    stop("Temporary FASTA file '", temp_fasta_path, "' ", size_info, ". Check permissions or available disk space.")
}

message("Starting Polyester simulation...")
message("Output will be written to: ", output_dir)

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Call the simulation function using the count matrix (all transcripts) AND the temporary FASTA file (all transcripts)
simulate_experiment_countmat(
    fasta = temp_fasta_path,        # Path to the temporary FASTA sequences (all used transcripts)
    readmat = count_matrix,         # The matrix of scaled mean counts (for all transcripts)
    outdir = output_dir,            # Directory to save FASTQ files
    readlen = READ_LENGTH,
    paired = PAIRED_END,
    fraglen = FRAGMENT_MEAN_LENGTH,
    fragsd = FRAGMENT_SD_LENGTH,
    # error_model = 'illumina5',    # Optional: specify an error model
    # num_threads = 4,              # Optional: use multiple cores if available
    seed = get(".Random.seed", envir = .GlobalEnv)[1] # Pass the current seed state
)

message("Polyester simulation finished.")
message("Simulated FASTQ files and simulation info are in: ", output_dir)

# The temporary file will be cleaned up automatically by the on.exit() call defined at the start

# --- End of Script ---
