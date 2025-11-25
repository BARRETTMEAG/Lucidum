
# THIS WORKS, pulled the core region needed from a supermatrix after using 
#Andrews Liger program.
############################################################
# Minimal Pipeline: Extract Core Region + View in msaR
# FIXED GAP HANDLING
############################################################

library(ape)
library(Biostrings)
library(msaR)

############################################################
# Step 1: Load FASTA
############################################################

fasta_file <- "smatrix.fasta"
aln <- read.dna(fasta_file, format = "fasta")

############################################################
# Step 2: Define TRUE GAP detection
############################################################

# DNAbin encodes:
# 0=A 1=C 2=G 3=T 4=N ... 240–255 = gap characters
is_gap <- function(x) {
  as.integer(x) >= 240   # TRUE gap
}

############################################################
# Step 3: Remove columns that are entirely gaps
############################################################

keep_cols <- apply(aln, 2, function(col) any(!is_gap(col)))
aln_no_gaps <- aln[, keep_cols, drop = FALSE]

############################################################
# Step 4: Identify core columns
# Core = columns with >=50% non-gaps
############################################################

num_non_gaps <- apply(aln_no_gaps, 2, function(col) sum(!is_gap(col)))
threshold <- 0.5 * nrow(aln_no_gaps)

core_cols <- num_non_gaps >= threshold

# Find longest continuous TRUE block
r <- rle(core_cols)

if (any(r$values)) {
  max_block <- which.max(r$lengths * r$values)
  end <- cumsum(r$lengths)[max_block]
  start <- end - r$lengths[max_block] + 1
  aln_core <- aln_no_gaps[, start:end, drop = FALSE]
} else {
  aln_core <- aln_no_gaps
}

cat("Core region columns:", start, "to", end, "\n")
cat("Core region length:", ncol(aln_core), "columns\n")

############################################################
# Step 5: Convert to DNAStringSet for msaR
############################################################

dnabin_to_DNAStringSet <- function(dna_matrix) {
  char_matrix <- as.character(dna_matrix)
  seqs <- apply(char_matrix, 1, paste, collapse = "")
  DNAStringSet(seqs, use.names = TRUE)
}

core_DNAStringSet <- dnabin_to_DNAStringSet(aln_core)

############################################################
# Step 6: View in msaR
############################################################

msaR(core_DNAStringSet)

############################################################
# Step 7: Save output
############################################################

writeXStringSet(core_DNAStringSet, "core_region.fasta")
cat("Core alignment saved to core_region.fasta\n")




# ##### This one was useful for trimming the core, but the one above is better#####
# ############################################################
# # Full Core Extraction + Diagnostics + msaR
# ############################################################
# 
# library(ape)
# library(Biostrings)
# library(msaR)
# 
# ############################################################
# # Step 1: Load FASTA
# ############################################################
# 
# fasta_file <- "smatrix.fasta"
# aln <- read.dna(fasta_file, format = "fasta")
# 
# ############################################################
# # Step 2: Define TRUE gap and N detection
# ############################################################
# 
# is_gap <- function(x) as.integer(x) >= 240   # real gaps
# is_N   <- function(x) as.integer(x) == 4     # N
# 
# ############################################################
# # Step 3: Remove all-gap columns
# ############################################################
# 
# keep_cols <- apply(aln, 2, function(col) any(!is_gap(col)))
# aln_ng <- aln[, keep_cols, drop = FALSE]
# 
# ############################################################
# # Step 4: Core region detection (>=50% non-gaps)
# ############################################################
# 
# num_non_gaps <- apply(aln_ng, 2, function(col) sum(!is_gap(col)))
# threshold <- 0.5 * nrow(aln_ng)
# 
# core_cols <- num_non_gaps >= threshold
# 
# r <- rle(core_cols)
# if (any(r$values)) {
#   max_block <- which.max(r$lengths * r$values)
#   end <- cumsum(r$lengths)[max_block]
#   start <- end - r$lengths[max_block] + 1
#   aln_core <- aln_ng[, start:end, drop = FALSE]
# } else {
#   aln_core <- aln_ng
# }
# 
# cat("Core region length:", ncol(aln_core), "columns\n")
# 
# ############################################################
# # Step 5: Remove taxa with any gap or N in core
# ############################################################
# 
# rows_to_keep <- apply(aln_core, 1, function(row) all(!is_gap(row) & !is_N(row)))
# 
# cat("Removed", sum(!rows_to_keep), "taxa due to internal gaps or N.\n")
# 
# aln_core_filtered <- aln_core[rows_to_keep, , drop = FALSE]
# 
# ############################################################
# # Step 6: DIAGNOSTICS — Raw DNAbin should be equal-length
# ############################################################
# 
# cat("\n=== DNAbin ROW LENGTH CHECK ===\n")
# print(apply(aln_core_filtered, 1, length))
# 
# # Convert each DNAbin byte to single characters (raw-safe)
# debug_strings <- apply(aln_core_filtered, 1, function(row) {
#   paste(rawToChar(row), collapse = "")
# })
# 
# cat("\n=== CHARACTER STRING LENGTHS (THESE MUST BE IDENTICAL) ===\n")
# print(nchar(debug_strings))
# 
# # Identify problematic rows
# bad_rows <- nchar(debug_strings) != ncol(aln_core_filtered)
# 
# if (any(bad_rows)) {
#   cat("\n!!! WARNING: THE FOLLOWING TAXA HAVE WRONG STRING LENGTHS !!!\n")
#   print(names(debug_strings)[bad_rows])
# }
# 
# ############################################################
# # Step 7: Guaranteed-safe DNAbin → DNAStringSet conversion
# ############################################################
# 
# dnabin_to_DNAStringSet <- function(dna_matrix) {
#   
#   # convert entire matrix to characters reliably
#   dna_chars <- matrix(
#     rawToChar(dna_matrix, multiple = TRUE),
#     nrow = nrow(dna_matrix),
#     byrow = TRUE
#   )
#   
#   # Diagnostics
#   cat("\n=== dna_chars ROW LENGTH CHECK ===\n")
#   print(apply(dna_chars, 1, nchar))
#   
#   seqs <- apply(dna_chars, 1, paste, collapse = "")
#   DNAStringSet(seqs, use.names = TRUE)
# }
# 
# core_DNAStringSet <- dnabin_to_DNAStringSet(aln_core_filtered)
# 
# ############################################################
# # Step 8: FINAL CHECK BEFORE msar
# ############################################################
# 
# cat("\n=== FINAL SEQUENCE LENGTH CHECK ===\n")
# print(nchar(as.character(core_DNAStringSet)))
# 
# if (length(unique(nchar(as.character(core_DNAStringSet)))) != 1) {
#   stop("ERROR: msaR cannot load because final sequences are not identical in length.")
# }
# 
# ############################################################
# # Step 9: View in msaR
# ############################################################
# 
# msaR(core_DNAStringSet)
# 
# ############################################################
# # Step 10: Save final FASTA
# ############################################################
# 
# writeXStringSet(core_DNAStringSet, "core_region.filtered.fasta")
# cat("Filtered core alignment saved to core_region.filtered.fasta\n")
