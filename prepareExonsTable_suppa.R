library(MatrixGenerics)
library(sparseMatrixStats)
library(DelayedMatrixStats)
###########################################
getwd()

dPSI_su = read.table("deltaSUPEx.tab")
dim(dPSI_su)
colnames(dPSI_su)
rownames(dPSI_su)<-dPSI_su$cleanSuppa.X2
dPSI_su$cleanSuppa.X1<- NULL
dPSI_su$cleanSuppa.X2<- NULL
dPSI_su<-dPSI_su*100
MIN <- rowMins(as.matrix(dPSI_su))
MAX <- rowMaxs(as.matrix(dPSI_su))
RANGE = MAX - MIN
dPSI_su_full = cbind(dPSI_su, RANGE = RANGE)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_su_full_chaging = dPSI_su_full[dPSI_su_full$RANGE > 5, ]
dim(dPSI_su_full_chaging) #1573 rows
#################################################

#so for the network I rename nodes later according to class table!!!
dPSI <- dPSI_su_full_chaging
dPSI$RANGE <- NULL
iis <- match(colnames(dPSI), class$Gene.Symbol)
iis

colnames(dPSI)[which(is.na(iis))]
# "AA1"       "AA2"       "AA3"       "AA4"       "AA5"       "AA6"       "AA7"      
# "AA8"       "AA9"       "C1orf55"   "C1orf55_b" "CWC22_b"   "HFM1_b"    "IKcon"    
#"PPIL1"     "PRPF8con"  "SF3B1con"  "SMU1con"   "SRPK2_b"   "SRRM4"     "SRRT"     
# "XAB2_b"
##########################################
##########################################
#single scaled
deltascaled <- scale(dPSI)# scaled by columns (KDs)
write.table(deltascaled, "ES_su_sscaled.tab",  sep = "\t")
####################################################
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
write.table(dsscaled, "ES_su_dscaled.tab",  sep = "\t")
####################################################


