#finalMatrixCOmparison
setwd("~/useful/")
library(data.table)
library(stringr)
INCFile<-""
OUTFile<-""
source("replaceNA.R")
replaceNA(INCFile, OUTFile)
#########################################
INCFile<-""
OUTFile1<-""
OUTFile2<-""
source("replaceNAI_corrected.R")
replaceNAI(INCFile, OUTFile1,OUTFile2)
#########################################
infile<-read.table("", 
                     header=T, sep="\t",na.strings = "NA" , row.names = 1)
head(infile)
dim(infile)
table(infile$COMPLEX)
nasC<-which(is.na(infile$COMPLEX))
write.table(infile[nasC,1:6], "", sep="\t", col.names = NA, quote = F)
rownames(infile)  <-paste(infile$EVENT, 1:nrow(infile), sep="_")
rownames(infile)[1:10]

