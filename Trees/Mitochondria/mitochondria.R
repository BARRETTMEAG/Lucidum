getwd()
library(ape)
library(phangorn)
library(msa)
library(seqinr)
library(msaR)

amy <- read.FASTA("Mitochondria/mito.fasta")

# Step 2: Check sequence count and unique names
length(amy)
length(unique(names(amy)))

# Step 3: Look at the MSA (if already aligned)
msaR(amy)

# Step 4: Remove duplicates
amy_unique <- amy[!duplicated(names(amy))]

# Step 5: Convert to phyDat for trimming
amy2 <- phyDat(amy_unique, type = "DNA")

# Step 6: Trim columns with >50% gaps
amy3 <- amy2[, colMeans(as.character(amy2) == "-") < 0.5]

# Step 7: Convert back to DNAbin and visualize
amy4 <- as.DNAbin(amy3)
msaR(amy4)

# Step 8: Save trimmed alignment to a new FASTA file
write.FASTA(amy4, "mito_trimmed.fasta")

amy1 <- read.FASTA("Mitochondria/Primates.fasta")

# Step 2: Check sequence count and unique names
length(amy1)
length(unique(names(amy1)))

# Step 3: Look at the MSA (if already aligned)
msaR(amy1)

# Step 4: Remove duplicates
amy_unique1 <- amy1[!duplicated(names(amy1))]

# Step 5: Convert to phyDat for trimming
amy2.1 <- phyDat(amy_unique1, type = "DNA")

# Step 6: Trim columns with >50% gaps
amy3.1 <- amy2.1[, colMeans(as.character(amy2.1) == "-") < 0.5]

# Step 7: Convert back to DNAbin and visualize
amy4.1 <- as.DNAbin(amy3.1)
msaR(amy4.1)

# Step 8: Save trimmed alignment to a new FASTA file
write.FASTA(amy4.1, "Primates_mito_trimmed.fasta")

#----------------------------------------
# Libraries
#----------------------------------------
library(seqinr)     # read/write FASTA
library(ape)        # DNAbin conversion
library(phangorn)   # phyDat & trimming
library(msaR)       # interactive MSA visualization

#----------------------------------------
# Step 1: Read trimmed FASTAs
#----------------------------------------
mito <- read.fasta("Mitochondria/mito_trimmed.fasta", as.string = TRUE)
primates <- read.fasta("Mitochondria/Primates_mito_trimmed.fasta", as.string = TRUE)

#----------------------------------------
# Step 2: Combine sequences
#----------------------------------------
combined <- c(mito, primates)
combined <- combined[!duplicated(names(combined))]  # remove duplicate names

#----------------------------------------
# Step 3: Convert to phyDat (already aligned)
#----------------------------------------
# Convert sequences to character vectors
combined_char <- lapply(combined, function(x) unlist(strsplit(x[[1]], split = "")))

# Convert to phyDat
combined_phy <- phyDat(combined_char, type = "DNA")

#----------------------------------------
# Step 4: Trim columns with >50% gaps
#----------------------------------------
gap_thresh <- 0.5
combined_trim <- combined_phy[, colMeans(as.character(combined_phy) == "-") < gap_thresh]

#----------------------------------------
# Step 5: Convert back to DNAbin for visualization
#----------------------------------------
combined_DNAbin <- as.DNAbin(combined_trim)

#----------------------------------------
# Step 6: Visualize MSA interactively
#----------------------------------------
msaR(combined_DNAbin)

#----------------------------------------
# Step 7: Optional - save trimmed combined alignment
#----------------------------------------
write.FASTA(combined_DNAbin, file.out = "combined_mito_trimmed.fasta")
