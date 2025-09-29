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