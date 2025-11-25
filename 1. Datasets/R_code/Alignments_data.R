######  MC1R  ###########
library(rentrez)
library(XML)
library(dplyr)
library(httr)
library(seqinr)

 

#attach search function for databases
entrez_dbs()
entrez_db_searchable(db="nuccore")
entrez_db_searchable(db="taxonomy")

 

#set_entrez_key()
#we want to do a search AND put it in an object
#Primate (Order) = txid9443
onesearch<-entrez_search(db="nuccore",term="Primate MC1R",
                         retmax=1)

 

twosearch<-entrez_search(db="nuccore", term="Primate[ORGN] MC1R[ALL]",
                         retmax=1)
twosearch$ids
#info on the search
onesearch$ids
onesearch$count
onesearch$retmax
onesearch$file

 

#we have a search result, but don't have specific info
sumsearch<-entrez_summary(db="nuccore", id=onesearch$ids)
sumsearch$slen

 

sumsearch2<-entrez_summary(db="nuccore", id=twosearch$ids)
sumsearch2$organisms

 

testsearch<-entrez_search(db="nuccore", term="txid9443[ORGN] AND MC1R[ALL]", 
                          retmax=1)
sumtest<-entrez_summary(db="nuccore", id=testsearch$ids) 
sumtest$organism
sumtest$accessionversion
sumsearch2$accessionversion

 

#getting sequence
output<-entrez_fetch(db="nuccore", id=testsearch$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="firstseq.fasta")

 

#we can get a lot of info from taxa
taxtest<-entrez_search(db="taxonomy", term="primate")
taxsum<-entrez_summary(db="taxonomy", id=taxtest$ids)
taxsum$taxid
taxout<-entrez_search(db="taxonomy", term="txid9443[SBTR]", retmax=2000)
rec<-entrez_summary(db="taxonomy", id=taxout$ids[3])
rec$scientificname

 

#for this to work, you need to have your taxonomy id output
#we need to make an empty table to put our data in
mc1rdat<-matrix(, nrow=(length(unique(taxout$ids))), ncol=5)
colnames(mc1rdat)<-c("taxid","speciesname","accnum","seqname","slen")
#What data do we want?
#taxid, speciesname, accnum, "seqname", slen
#i<-42
mitosequence<-character()
for(i in 1:length(unique(taxout$ids))){
  print(i)
  mc1rdat[i,1]<-taxout$ids[i]
  record<-entrez_summary(db="taxonomy", id=taxout$ids[i])
  mc1rdat[i,2]<-record$scientificname
  seqout<-entrez_search(db="nuccore", term=paste("txid", taxout$ids[i],"[ORGN] 
                        AND MC1R[ALL] AND 1:10000 [SLEN]",
                        sep=""), retmax=1)
  if(length(seqout$ids)<1){
    mc1rdat[i,3]<-NA
    mc1rdat[i,4]<-NA
    mc1rdat[i,5]<-NA
  }
  else{
    sumout<-entrez_summary(db="nuccore", id=seqout$ids)
    mc1rdat[i,3]<-sumout$accessionversion
    mc1rdat[i,4]<-sumout$title
    mc1rdat[i,5]<-sumout$slen
    seq.dat<-entrez_fetch(db="nuccore", id=sumout$accessionversion, rettype="fasta")
    seq.dat<-sub(">([^\n]*)", paste0(">", sumout$organism), seq.dat)
    seq.dat<-gsub(" ", "_", seq.dat)
    mitosequence<-c(mitosequence, seq.dat)
  }
}
write.csv(mc1rdat, "MC1R.csv")
write(mitosequence, file="MC1R_out.fasta")

####### RAG2 #########
library(rentrez)
library(XML)
library(dplyr)
library(httr)
library(seqinr)



#attach search function for databases
entrez_dbs()
entrez_db_searchable(db="nuccore")
entrez_db_searchable(db="taxonomy")



#set_entrez_key()
#we want to do a search AND put it in an object
#Primate (Order) = txid9443
onesearch<-entrez_search(db="nuccore",term="Primate RAG2",
                         retmax=1)



twosearch<-entrez_search(db="nuccore", term="Primate[ORGN] RAG2[ALL]",
                         retmax=1)
twosearch$ids
#info on the search
onesearch$ids
onesearch$count
onesearch$retmax
onesearch$file



#we have a search result, but don't have specific info
sumsearch<-entrez_summary(db="nuccore", id=onesearch$ids)
sumsearch$slen



sumsearch2<-entrez_summary(db="nuccore", id=twosearch$ids)
sumsearch2$organisms



testsearch<-entrez_search(db="nuccore", term="txid9443[ORGN] AND RAG2[ALL]", 
                          retmax=1)
sumtest<-entrez_summary(db="nuccore", id=testsearch$ids) 
sumtest$organism
sumtest$accessionversion
sumsearch2$accessionversion



#getting sequence
output<-entrez_fetch(db="nuccore", id=testsearch$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="firstseq1.fasta")



#we can get a lot of info from taxa
taxtest<-entrez_search(db="taxonomy", term="primate")
taxsum<-entrez_summary(db="taxonomy", id=taxtest$ids)
taxsum$taxid
taxout<-entrez_search(db="taxonomy", term="txid9443[SBTR]", retmax=2000)
rec<-entrez_summary(db="taxonomy", id=taxout$ids[3])
rec$scientificname



#for this to work, you need to have your taxonomy id output
#we need to make an empty table to put our data in
mc1rdat<-matrix(, nrow=(length(unique(taxout$ids))), ncol=5)
colnames(mc1rdat)<-c("taxid","speciesname","accnum","seqname","slen")
#What data do we want?
#taxid, speciesname, accnum, "seqname", slen
#i<-42
mitosequence1<-character()
for(i in 1:length(unique(taxout$ids))){
  print(i)
  mc1rdat[i,1]<-taxout$ids[i]
  record<-entrez_summary(db="taxonomy", id=taxout$ids[i])
  mc1rdat[i,2]<-record$scientificname
  seqout<-entrez_search(db="nuccore", term=paste("txid", taxout$ids[i],"[ORGN] 
                        AND RAG2[ALL] AND 1:10000 [SLEN]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    mc1rdat[i,3]<-NA
    mc1rdat[i,4]<-NA
    mc1rdat[i,5]<-NA
  }
  else{
    sumout<-entrez_summary(db="nuccore", id=seqout$ids)
    mc1rdat[i,3]<-sumout$accessionversion
    mc1rdat[i,4]<-sumout$title
    mc1rdat[i,5]<-sumout$slen
    seq.dat<-entrez_fetch(db="nuccore", id=sumout$accessionversion, rettype="fasta")
    seq.dat<-sub(">([^\n]*)", paste0(">", sumout$organism), seq.dat)
    seq.dat<-gsub(" ", "_", seq.dat)
    mitosequence1<-c(mitosequence1, seq.dat)
  }
}
write.csv(mc1rdat, "RAG2.csv")
write(mitosequence1, file="RAG2_out.fasta")



############################################################
## Process and Trim Individual Nucleotide Files (and save as .fasta)
############################################################

# --- RAG2 ---
amy_dna <- readDNAStringSet("RAG2_out.fasta")
amy_dna <- amy_dna[!duplicated(names(amy_dna))] 
names(amy_dna) <- gsub(" ", "_", names(amy_dna)) 
aligned <- msa(amy_dna, method = "ClustalW")
aligned_dna <- as(aligned, "DNAStringSet")
aligned_matrix <- as.matrix(aligned_dna)
amy_phy <- phyDat(aligned_matrix, type = "DNA") 
gap_mask <- colMeans(aligned_matrix == "-") < 0.5
amy_trim <- amy_phy[, gap_mask]
amy_trim_matrix <- as.character(amy_trim)
amy_trim_dna <- DNAStringSet(sapply(1:nrow(amy_trim_matrix), function(i) paste(amy_trim_matrix[i, ], collapse = "")))
names(amy_trim_dna) <- rownames(amy_trim_matrix)
writeXStringSet(amy_trim_dna, "RAG2_trimmed.fasta") 

# --- MC1R ---
amy_dna2 <- readDNAStringSet("MC1R_out.fasta")
amy_dna2 <- amy_dna2[!duplicated(names(amy_dna2))] 
names(amy_dna2) <- gsub(" ", "_", names(amy_dna2)) 
aligned2 <- msa(amy_dna2, method = "ClustalW")
aligned_dna2 <- as(aligned2, "DNAStringSet")
aligned_matrix2 <- as.matrix(aligned_dna2)
amy_phy2 <- phyDat(aligned_matrix2, type = "DNA")
gap_mask2 <- colMeans(aligned_matrix2 == "-") < 0.5
amy_trim2 <- amy_phy2[, gap_mask2]
amy_trim_matrix2 <- as.character(amy_trim2)
amy_trim_dna2 <- DNAStringSet(sapply(1:nrow(amy_trim_matrix2), function(i) paste(amy_trim_matrix2[i, ], collapse = "")))
names(amy_trim_dna2) <- rownames(amy_trim_matrix2)
writeXStringSet(amy_trim_dna2, "MC1R_trimmed.fasta") 

# --- IRBP ---
amy_dna3 <- readDNAStringSet("IRBP_trimmed2.fasta")
amy_dna3 <- amy_dna3[!duplicated(names(amy_dna3))] 
names(amy_dna3) <- gsub(" ", "_", names(amy_dna3)) 
aligned3 <- msa(amy_dna3, method = "ClustalW")
aligned_dna3 <- as(aligned3, "DNAStringSet")
aligned_matrix3 <- as.matrix(aligned_dna3)
amy_phy3 <- phyDat(aligned_matrix3, type = "DNA")
gap_mask3 <- colMeans(aligned_matrix3 == "-") < 0.5
amy_trim3 <- amy_phy3[, gap_mask3]
amy_trim_matrix3 <- as.character(amy_trim3)
amy_trim_dna3 <- DNAStringSet(sapply(1:nrow(amy_trim_matrix3), function(i) paste(amy_trim_matrix3[i, ], collapse = "")))
names(amy_trim_dna3) <- rownames(amy_trim_matrix3)
writeXStringSet(amy_trim_dna3, "IRBP_trimmed.fasta") 

# --- mtDNA ---
amy_dna4 <- readDNAStringSet("mtDNA_trimmed.fasta")
amy_dna4 <- amy_dna4[!duplicated(names(amy_dna4))] 
names(amy_dna4) <- gsub(" ", "_", names(amy_dna4)) 
aligned4 <- msa(amy_dna4, method = "ClustalW")
aligned_dna4 <- as(aligned4, "DNAStringSet")
aligned_matrix4 <- as.matrix(aligned_dna4)
amy_phy4 <- phyDat(aligned_matrix4, type = "DNA")
gap_mask4 <- colMeans(aligned_matrix4 == "-") < 0.5
amy_trim4 <- amy_phy4[, gap_mask4]
amy_trim_matrix4 <- as.character(amy_trim4)
amy_trim_dna4 <- DNAStringSet(sapply(1:nrow(amy_trim_matrix4), function(i) paste(amy_trim_matrix4[i, ], collapse = "")))
names(amy_trim_dna4) <- rownames(amy_trim_matrix4)
writeXStringSet(amy_trim_dna4, "mtDNA_trimmed1.fasta") 

# --- prim_GPD6 ---
amy_dna5 <- readDNAStringSet("prim_GPD6_Trim.fasta")
amy_dna5 <- amy_dna5[!duplicated(names(amy_dna5))] 
names(amy_dna5) <- gsub(" ", "_", names(amy_dna5)) 
aligned5 <- msa(amy_dna5, method = "ClustalW")
aligned_dna5 <- as(aligned5, "DNAStringSet")
aligned_matrix5 <- as.matrix(aligned_dna5)
amy_phy5 <- phyDat(aligned_matrix5, type = "DNA")
gap_mask5 <- colMeans(aligned_matrix5 == "-") < 0.5
amy_trim5 <- amy_phy5[, gap_mask5]
amy_trim_matrix5 <- as.character(amy_trim5)
amy_trim_dna5 <- DNAStringSet(sapply(1:nrow(amy_trim_matrix5), function(i) paste(amy_trim_matrix5[i, ], collapse = "")))
names(amy_trim_dna5) <- rownames(amy_trim_matrix5)
writeXStringSet(amy_trim_dna5, "prim_GPD6_Trim1.fasta") 


#Then used Liger to create supermatrix. View Core_region.R for more notes