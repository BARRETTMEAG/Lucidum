####We are trying to get the mitochondrion of members of the Hominoidae family####
#libraries/packages needed
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
#Hominoidea (Superfamily) = txid314295
#Hominidae (Family) = txid9604
onesearch<-entrez_search(db="nuccore",term="Hominidae mitochondrion",
                         retmax=1)

twosearch<-entrez_search(db="nuccore", term="Hominidae[ORGN] mitochondrion[ALL]",
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

testsearch<-entrez_search(db="nuccore", term="txid9604[ORGN] and MITOCHONDRION[ALL]", 
                          retmax=1) #Why is it saying there are no ids for Hominidae or Hominoidea?
sumtest<-entrez_summary(db="nuccore", id=testsearch$ids) 
sumtest$organism
sumtest$accessionversion
sumsearch2$accessionversion

#getting sequence
output<-entrez_fetch(db="nuccore", id=testsearch$ids, rettype="fasta")
cat(strwrap(substr(output, 1, 500)), sep="\n")
write(output, file="firstseq.fasta")

#we can get a lot of info from taxa
taxtest<-entrez_search(db="taxonomy", term="Hominidae")
taxsum<-entrez_summary(db="taxonomy", id=taxtest$ids)
taxsum$taxid
taxout<-entrez_search(db="taxonomy", term="txid9604[SBTR]", retmax=2000)
rec<-entrez_summary(db="taxonomy", id=taxout$ids[3])
rec$scientificname

#for this to work, you need to have your taxonomy id output
#we need to make an empty table to put our data in
homidat<-matrix(, nrow=(length(unique(taxout$ids))), ncol=5)
colnames(homidat)<-c("taxid","speciesname","accnum","seqname","slen")
#What data do we want?
#taxid, speciesname, accnum, "seqname", slen
#i<-42
mitosequence<-character()
for(i in 1:length(unique(taxout$ids))){
  print(i)
  homidat[i,1]<-taxout$ids[i]
  record<-entrez_summary(db="taxonomy", id=taxout$ids[i])
  homidat[i,2]<-record$scientificname  #Why isn't this line working?
  seqout<-entrez_search(db="nuccore", term=paste("txid", taxout$ids[i],"[ORGN] AND MITOCHONDRION[ALL]", sep=""),
                        retmax=1)
  if(length(seqout$ids)<1){
    homidat[i,3]<-NA
    homidat[i,4]<-NA
    homidat[i,5]<-NA
  }
  else{
    sumout<-entrez_summary(db="nuccore", id=seqout$ids)
    homidat[i,3]<-sumout$accessionversion
    homidat[i,4]<-sumout$title
    homidat[i,5]<-sumout$slen
    seq.dat<-entrez_fetch(db="nuccore", id=sumout$accessionversion, rettype="fasta")
    seq.dat<-sub(">([^\n]*)", paste0(">", sumout$organism), seq.dat) #to simplify sequence name
    seq.dat<-gsub(" ", "_", seq.dat)
    mitosequence<-c(mitosequence, seq.dat)
  }
}
write(mitosequence, file="mitochondrial_out.fasta")
