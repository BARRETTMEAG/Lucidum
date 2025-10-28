#load packages
#install.packages("evobiR")
#install.packages("msaR")
library(evobiR)
library(ape)
library(msaR)
library(phangorn)

homi<-read.FASTA("Homi_mito_alignment.fasta")#put in fasta alignment
length(homi)

#look at MSA
msaR(homi)
length(unique(names(homi)))
#remove duplicates
homi_unique<-homi[!duplicated(names(homi))]
#trimming
homi2<-phyDat(homi_unique, type="DNA")
#this trims data
homi3<-homi2[, colMeans(as.character(homi2)== "-")<0.5]
homi4<-as.DNAbin(homi3)
msaR(homi4)
write.FASTA(homi4,"Homi_mito_trimmed.fasta")


#making a supermatrix that will run all fasta in working directory
SuperMatrix(missing="-", prefix = "Homi_sm", save=T)

#new package
#install.packages("Rogue")
library(Rogue)
boots<-read.tree("Homi_mito_trimmed.fasta.ufboot")
#plot (boots[1])
testing<-TipInstability(boots)
hist(testing)
contree<-read.tree("Homi_mito_trimmed.fasta.contree")
plot(contree, cex=0.4)
plot(contree, tip.col=ColByStability(boots), cex=0.4)
names(testing[testing>0.04]) #change the X.XX as needed
