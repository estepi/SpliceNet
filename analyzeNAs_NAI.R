#finalMatrixCOmparison
library(data.table)
library(stringr)
#########################################
#replace new NAs and NAIs at once
INCFile<-""
OUTFile1<-""
OUTFile2<-""
source("replaceNAI.R")
replaceNAI(INCFile, OUTFile1,OUTFile2)
#########################################
# Only replace new NAs
INCFile<-""
OUTFile<-""
source("replaceNA.R")
replaceNA(INCFile, OUTFile)

