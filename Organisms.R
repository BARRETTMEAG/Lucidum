# load files
library(rentrez)
library(XML)
library(dplyr)
library(httr)
library(seqinr)


# 1st sequence Pongo pygmaeus in nuccore
firstsearch1 <-entrez_search(db="nuccore", term= "Pongo pygmaeus", 
                             retmax=1)
firstsearch1$ids
firstsearch1$count
firstsearch1$retmax
firstsearch1$file

# Searching all pongo pygmaeus mitochondrion 
secsearch2 <-entrez_search(db = "nuccore", term = "Pongo pygmaeus[ORGN] AND MITOCHONDRION[ALL]", 
                           retmax =1)

sumsearch2 <-entrez_summary(db = "nuccore", id=firstsearch1$ids)
sumsearch1$slen

secsearch2$ids
sumsearch2.1 <- entrez_summary(db= "nuccore", id=secsearch2$ids)
sumsearch2.1$organism
testsearch1 <-entrez_search(db = "nuccore", term="txid9600[ORGN] AND MITOCHONDRION[ALL]", retmax=1)


#we have a search result, but don't have specific info
sumsearch1 <-entrez_summary(db = "nuccore", id=firstsearch1$ids)
sumsearch1$slen

secsearch2$ids
sumsearch2.1 <- entrez_summary(db= "nuccore", id=secsearch2$ids)
sumsearch2.1$organism

testsearch1 <-entrez_search(db = "nuccore", term="txid9600[ORGN] AND MITOCHONDRION[ALL]", retmax=1)
testsearch1

sumtest1 <-entrez_summary(db="nuccore", id=testsearch1$ids)
sumtest1$organism
sumtest1$accessionversion


# writing out the file in fasta
output1 <- entrez_fetch(db="nuccore",id=testsearch1$ids, rettype="fasta")
cat(strwrap(substr(output1, 1,500)),sept="\n")
write(output1, file="firstseq1.fasta")

taxtest1<-entrez_search(db="taxonomy", term="Pongo")
taxout1<-entrez_search(db="taxonomy", term="txid9600[SBTR]", retmax= 2000)
#taxout2<-entrez_search(db="taxonomy", term="txid2753606[SBTR] AND MITOCHONDRION[ALL]", retmax= 2000)
rec1 <-entrez_summary(db="taxonomy", id=taxout1$ids[1])
rec1$scientificname
rec1$commonname

mdat<-matrix(, nrow=(length(unique(taxout1$ids))), ncol =5)
colnames(mdat)<- c("txid", "speciesname","accnum", "seqname", "slen")

#think about what we want
#txid, speciesname, accum, "seqname", slen (sequence name)
i<-1
m_sequence<-character()
for(i in 1:length(unique(taxout1$ids))){
  print(i)
  mdat[i,1]<-taxout1$ids[i]
  record_0<-entrez_summary(db="taxonomy", id=taxout1$ids[i])
  mdat[i,2]<-record_0$scientificname
  se_out<-entrez_search(db="nuccore", term=paste("txid", taxout1$ids[i], "[ORGN] AND MITOCHONDRION[ALL]", 
                                                  sep = ""), retmax = 1) 
  if(length(se_out$ids)<1){
    mdat[i,3]<-NA
    mdat[i,4]<-NA
    mdat[i,5]<-NA
  }
  else{
    s_out<-entrez_summary(db="nuccore", id=se_out$ids)
    mdat[i,3]<-s_out$accessionversion
    mdat[i,4]<-s_out$title
    mdat[i,5]<-s_out$slen
    se_1.2<-entrez_fetch(db="nuccore", id=s_out$accessionverstion,rettype="fasta")
    se_dat<-sub(">([^\n]*", paste0(">", s_out$organism, se_dat))
    se_dat<-gsub(" ", "_")
    mi_sequence<-c(m_sequence, se_dat)
  }
}

# 26 Sept 25
write(mi_sequence, file = "p_mitochondrial_out.fasta")

# Cercopithecidae
newsearch <-entrez_search(db="nuccore", term= "Cercopithecidae", 
                             retmax=1)
newsearch$ids
newsearch$count
newsearch$retmax
newsearch$file

# Searching all Cercopithecidae mitochondrion 
newsearch2 <-entrez_search(db = "nuccore", term = "Cercopithecidae[ORGN] AND MITOCHONDRION[ALL]", 
                           retmax =1)

sum_search <-entrez_summary(db = "nuccore", id=newsearch$ids)
sum_search$slen

newsearch2$ids
sum_search2 <- entrez_summary(db= "nuccore", id=newsearch2$ids)
sum_search2$organism
test_search1 <-entrez_search(db = "nuccore", term="txid9527[ORGN] AND MITOCHONDRION[ALL]", retmax=1)

#we have a search result, but don't have specific info
sum_search1 <-entrez_summary(db = "nuccore", id=newsearch$ids)
sum_search1$slen

newsearch2$ids
sum_search2.1 <- entrez_summary(db= "nuccore", id=newsearch2$ids)
sum_search2.1$organism

test_search1 <-entrez_search(db = "nuccore", term="txid9527[ORGN] AND MITOCHONDRION[ALL]", retmax=1)
test_search1

sum_test1 <-entrez_summary(db="nuccore", id=test_search1$ids)
sum_test1$organism
sum_test1$accessionversion


# writing out the file in fasta
output2 <- entrez_fetch(db="nuccore",id=test_search1$ids, rettype="fasta")
cat(strwrap(substr(output2, 1,500)),sept="\n")
write(output2, file="secseq1.fasta")

tax_test1<-entrez_search(db="taxonomy", term="Cercopithecidae")
tax_out1<-entrez_search(db="taxonomy", term="txid9527[SBTR]", retmax= 2000)
rec_1 <-entrez_summary(db="taxonomy", id=tax_out1$ids[1])
rec_1$scientificname
rec_1$commonname

dat<-matrix(, nrow=(length(unique(tax_out1$ids))), ncol =5)
colnames(dat)<- c("txid", "speciesname","accnum", "seqname", "slen")

#think about what we want
#txid, speciesname, accum, "seqname", slen (sequence name)
i<-1
mito_sequence<-character()
for(i in 1:length(unique(tax_out1$ids))){
  print(i)
  dat[i,1]<-tax_out1$ids[i]
  record_1<-entrez_summary(db="taxonomy", id=tax_out1$ids[i])
  dat[i,2]<-record_1$scientificname
  seq_out<-entrez_search(db="nuccore", term=paste("txid", tax_out1.2$ids[i], "[ORGN] AND MITOCHONDRION[ALL]", 
                                                  sep = ""), retmax = 1) 
  if(length(seq_out$ids)<1){
    dat[i,3]<-NA
    dat[i,4]<-NA
    dat[i,5]<-NA
  }
  else{
    sum_out<-entrez_summary(db="nuccore", id=seq_out$ids)
    dat[i,3]<-sumout$accessionversion
    dat[i,4]<-sumout$title
    dat[i,5]<-sumout$slen
    seq_1<-entrez_fetch(db="nuccore", id=sum_out$accessionverstion,rettype="fasta")
    seq_1.dat<-sub(">([^\n]*", paste0(">", sum_out$organism, seq.dat))
    seq_1.dat<-gsub(" ", "_")
    mito_sequence<-c(mito_sequence, seq_dat)
  }
}

# 26 Sept 25
write(mito_sequence, file = "c_mitochondrial_out.fasta")

# Petromyzon marinus
new_search1 <-entrez_search(db="nuccore", term= "Petromyzon marinus", 
                          retmax=1)
new_search1$ids
new_search1$count
new_search1$retmax
new_search1$file

# Searching all Petromyzon marinus mitochondrion 
new_search2 <-entrez_search(db = "nuccore", term = "Petromyzon marinus[ORGN] AND MITOCHONDRION[ALL]", 
                           retmax =1)

sum_search1 <-entrez_summary(db = "nuccore", id=new_search1$ids)
sum_search1$slen

new_search2$ids
sum_search2.0 <- entrez_summary(db= "nuccore", id=new_search2$ids)
sum_search2.0$organism
test_search1.1 <-entrez_search(db = "nuccore", term="txid7757[ORGN] AND MITOCHONDRION[ALL]", retmax=1)

#we have a search result, but don't have specific info
sum_search1.2 <-entrez_summary(db = "nuccore", id=new_search1$ids)
sum_search1.2$slen

new_search2$ids
sum_search2.2 <- entrez_summary(db= "nuccore", id=newsearch2$ids)
sum_search2.1$organism

test_search1.1 <-entrez_search(db = "nuccore", term="txid7757[ORGN] AND MITOCHONDRION[ALL]", retmax=1)
test_search1.1

sum_test1.0 <-entrez_summary(db="nuccore", id=test_search1.1$ids)
sum_test1.0$organism
sum_test1.0$accessionversion


# writing out the file in fasta
output3 <- entrez_fetch(db="nuccore",id=test_search1.1$ids, rettype="fasta")
cat(strwrap(substr(output3, 1,500)),sept="\n")
write(output3, file="thirdseq1.fasta")

tax_test1.2<-entrez_search(db="taxonomy", term="Petromyzon marinus")
tax_out1.2<-entrez_search(db="taxonomy", term="txid7757[SBTR]", retmax= 2000)
rec_1.1 <-entrez_summary(db="taxonomy", id=tax_out1.2$ids[1])
rec_1.1$scientificname
rec_1.1$commonname

summary(tax_out1.2)
tax_out1.2$ids

qdat<-matrix(, nrow=(length(unique(tax_out1.2$ids))), ncol =5)
colnames(dat)<- c("txid", "speciesname","accnum", "seqname", "slen")

#think about what we want
#txid, speciesname, accum, "seqname", slen (sequence name)
i<-1
mit_sequence<-character()
for(i in 1:length(unique(tax_out1.2$ids))){
  print(i)
  qdat[i,1]<-tax_out1.2$ids[i]
  record_1<-entrez_summary(db="taxonomy", id=tax_out1.2$ids[i])
  qdat[i,2]<-record_1$scientificname
  seq_out<-entrez_search(db="nuccore", term=paste("txid", tax_out1.2$ids[i], "[ORGN] AND MITOCHONDRION[ALL]", 
                                                 sep = ""), retmax = 1) 
  if(length(seq_out$ids)<1){
    dat[i,3]<-NA
    dat[i,4]<-NA
    dat[i,5]<-NA
  }
  else{
    sum_out<-entrez_summary(db="nuccore", id=seq_out$ids)
    dat[i,3]<-sumout$accessionversion
    dat[i,4]<-sumout$title
    dat[i,5]<-sumout$slen
    seq_1.2<-entrez_fetch(db="nuccore", id=sum_out$accessionverstion,rettype="fasta")
    seq_1.2.dat<-sub(">([^\n]*", paste0(">", sum_out$organism, seq.dat))
    seq_1.2dat<-gsub(" ", "_")
    mit_sequence<-c(mit_sequence, seq_1.2.dat)
  }
}

# 26 Sept 25
write(mit_sequence, file = "pm_mitochondrial_out.fasta")

