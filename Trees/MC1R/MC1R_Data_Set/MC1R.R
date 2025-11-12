####Doesn't Work####
# amy<-read.FASTA("MC1R_refseq_transcript.fasta")
# length(amy)
# length(unique(names(amy)))
# #look at the MSA
# msaR(amy)
# length(unique(names(amy)))
# #trim out duplicates
# any_unique <- amy[!duplicated(names(amy))]
# #trimming
# amy2<-phyDat(amyPunique, type="DNA")
# #this is where we trim out data, 80%
# amy3<-amy2[, colMeans(as.character(amy2)== "-")< 0.5]
# amy4 <-as.DNAbin(amy3)
# msaR(amy4)
# write.FASTA(amy4, "MC1R.fasta")
# 
# MC1R1<-read.FASTA("MC1R.fasta")
# length(MC1R1)
# msaR(MC1R1)
# MC1R1_unique<-MC1R1[!duplicated(names(MC1R1))]
# MC1R2<-phyDat(MC1R1_unique, type="DNA")
# MC1R3<-MC1R2[,colMeans(as.character(MC1R2)== "-")<0.5]
# MC1R4<-as.DNAbin(MC1R3)
# msaR(MC1R4)


# ------------------------------
# 1️⃣ Load required libraries
# ------------------------------
library(seqinr)       # For reading/writing FASTA
library(Biostrings)   # For AAStringSet and alignment handling
library(msa)          # For multiple sequence alignment
library(msaR)         # For visualizing alignments
library(phangorn)     # For phyDat conversion and trimming

# ------------------------------
# 2️⃣ Read protein FASTA
# ------------------------------
# Using Biostrings to directly get AA sequences
amy_aa <- readAAStringSet("MC1R.fasta")

# Remove duplicate sequence names
amy_aa <- amy_aa[!duplicated(names(amy_aa))]

# ------------------------------
# 3️⃣ Optional: visualize raw sequences
# ------------------------------
# Note: msaR requires aligned sequences
# So we'll align first before using msaR

# ------------------------------
# 4️⃣ Align sequences if not already aligned
# ------------------------------
aligned <- msa(amy_aa, method = "ClustalW")  # or method = "Muscle"

# Convert to AAStringSet for visualization
aligned_aa <- as(aligned, "AAStringSet")

# ------------------------------
# 5️⃣ Visualize alignment (optional)
# ------------------------------
msaR(aligned_aa)

# ------------------------------
# 6️⃣ Convert alignment to phangorn phyDat (AA)
# ------------------------------
aligned_matrix <- as.matrix(aligned_aa)
amy_phy <- phyDat(aligned_matrix, type = "AA")

# ------------------------------
# 7️⃣ Trim columns with >50% gaps
# ------------------------------
gap_mask <- colMeans(aligned_matrix == "-") < 0.5
amy_trim <- amy_phy[, gap_mask]

# ------------------------------
# 8️⃣ Convert trimmed phyDat back to character matrix
# ------------------------------
amy_trim_matrix <- as.character(amy_trim)

# ------------------------------
# 9️⃣ Write trimmed alignment to FASTA
# ------------------------------
write.fasta(
  sequences = lapply(1:nrow(amy_trim_matrix), function(i) amy_trim_matrix[i, ]),
  names = rownames(amy_trim_matrix),
  file.out = "MC1R_protein_trimmed.fasta"
)

# ------------------------------
# ✅ Done: trimmed protein alignment saved
# ------------------------------
# You can now visualize again if needed
amy_trim_aa <- AAStringSet(sapply(1:nrow(amy_trim_matrix), 
                                  function(i) paste(amy_trim_matrix[i, ], collapse = "")))
names(amy_trim_aa) <- rownames(amy_trim_matrix)
msaR(amy_trim_aa)

