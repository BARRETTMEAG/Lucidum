###############################################################
# FULL PHYLOGENETIC PIPELINE IN R
# ML + BOOTSTRAP + BAYESIAN BOOTSTRAP + ROOTING WITH Cynocephalus_volans
###############################################################

### 1. INSTALL & LOAD PACKAGES ################################
packages <- c("ape", "phangorn", "DECIPHER", "seqinr")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p, dependencies = TRUE)
}

library(ape)
library(phangorn)
library(DECIPHER)
library(seqinr)

###############################################################
### 2. READ FASTA FILE ########################################
###############################################################

fasta_file <- "new_core_region.fasta"

cat("Reading FASTA file...\n")
raw_sequences <- readDNAStringSet(fasta_file)
print(raw_sequences)

###############################################################
### 3. REMOVE GAPS & ALIGN SEQUENCES ##########################
###############################################################

cat("Removing gaps...\n")
raw_nogaps <- DNAStringSet(gsub("-", "", as.character(raw_sequences)))
names(raw_nogaps) <- names(raw_sequences)

cat("Aligning sequences with DECIPHER...\n")
alignment <- AlignSeqs(raw_nogaps, processors = NULL)

aligned_fasta <- "aligned_core_region.fasta"
writeXStringSet(DNAStringSet(alignment), aligned_fasta)
cat("Alignment saved to:", aligned_fasta, "\n")

###############################################################
### 4. ML PHYLOGENY (GTR+G+I) #################################
###############################################################

cat("Preparing data for ML analysis...\n")
phydat <- phyDat(as.matrix(alignment), type = "DNA")

dm <- dist.ml(phydat)
nj_tree <- NJ(dm)

cat("Optimizing ML tree with GTR+G+I...\n")
fit <- pml(nj_tree, phydat)

fit_GTR <- optim.pml(
  fit,
  model = "GTR",
  optNni = TRUE,
  optGamma = TRUE,
  optInv = TRUE,
  optBf = TRUE,
  optQ = TRUE
)

###############################################################
### 5. ML BOOTSTRAP ###########################################
###############################################################

cat("Running 100 ML bootstrap replicates...\n")
bs <- bootstrap.pml(
  fit_GTR,
  bs = 100,
  optNni = TRUE,
  multicore = FALSE
)

ML_tree <- plotBS(fit_GTR$tree, bs, p = 50)

###############################################################
### 6. ROOT ML TREE ###########################################
###############################################################

outgroup <- "Cynocephalus_volans"

cat("Rooting ML tree with:", outgroup, "\n")
rooted_ML <- root(ML_tree, outgroup = outgroup, resolve.root = TRUE)

write.tree(rooted_ML, "ML_tree_rooted.nwk")

pdf("ML_tree_rooted.pdf", width = 12, height = 12)
plot(rooted_ML, cex = 0.5)
title("Maximum Likelihood Tree (rooted on Cynocephalus_volans)")
dev.off()

###############################################################
### 7. BAYESIAN BOOTSTRAP (phangorn) ##########################
###############################################################

cat("Running Bayesian bootstrap (phangorn)...\n")

bb_trees <- bootstrap.phyDat(
  phydat,
  FUN = function(x) {
    tree <- optim.pml(
      pml(nj_tree, x),
      model = "GTR",
      optNni = TRUE,
      optGamma = TRUE,
      optInv = TRUE
    )$tree
    return(tree)
  },
  bs = 100
)

###############################################################
### 8. CONSENSUS TREE #########################################
###############################################################

cat("Computing Bayesian bootstrap consensus tree...\n")

# Majority-rule consensus using ape
bb_consensus <- ape::consensus(bb_trees)

###############################################################
### 9. ROOT BAYESIAN CONSENSUS TREE ###########################
###############################################################

cat("Rooting Bayesian bootstrap consensus tree...\n")

rooted_bb <- root(bb_consensus, outgroup = outgroup, resolve.root = TRUE)

write.tree(rooted_bb, "BayesBootstrap_tree_rooted.nwk")

###############################################################
### 10. PLOT BAYESIAN TREE ####################################
###############################################################

pdf("BayesBootstrap_tree_rooted.pdf", width = 12, height = 12)
plot(rooted_bb, cex = 0.5)
title("Bayesian Bootstrap Consensus Tree\n(rooted on Cynocephalus_volans)")
dev.off()
