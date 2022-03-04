#Estefania March 20222
#This script counts how many NAs each PSI value has.
#It PSi value has more than 3 NAs it will be replace for NAnew
# For IR, if pVal of READ UNBLANCE is significant (<0.05) PIR value will be replaced by NAI
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

