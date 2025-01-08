library(MatrixGenerics)
library(sparseMatrixStats)
library(DelayedMatrixStats)
###########################################
getwd()

dPSI_vt = read.table("deltaVTEx.tab")
dim(dPSI_vt)
colnames(dPSI_vt)
rownames(dPSI_vt)<-dPSI_vt$EVENT
dPSI_vt$EVENT<- NULL
dPSI_vt[1:10, 1:10]

MIN <- rowMins(as.matrix(dPSI_vt))
MAX <- rowMaxs(as.matrix(dPSI_vt))
RANGE = MAX - MIN
dPSI_vt_full = cbind(dPSI_vt, RANGE = RANGE)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_vt_full_chaging = dPSI_vt_full[dPSI_vt_full$RANGE > 5, ]
dim(dPSI_vt_full_chaging) #1348 rows
#################################################
colnames(dPSI_vt_full_chaging)[colnames(dPSI_vt_full_chaging) == "LENG1_b"]<-
  "LENG1"

colnames(dPSI_vt_full_chaging)[colnames(dPSI_vt_full_chaging) == "RBM17con"]<-
  "RBM17"

colnames(dPSI_vt_full_chaging)[colnames(dPSI_vt_full_chaging) == "HFM1_b"] <-  "HFM1"


colnames(dPSI_vt_full_chaging)[colnames(dPSI_vt_full_chaging) == "CCDC12_b"] <-
  "CCDC12"

colnames(dPSI_vt_full_chaging)[colnames(dPSI_vt_full_chaging) == "CDC5L_b"] <-
  "CDC5L"

dim(dPSI_vt_full_chaging)
#original table is names according to Gene.Symbol

class <-
  read.delim(
    "/Users/estepi/Documents/CRG/SpliceNetData/class_colors_2020.txt",
header = T)

#so for the network I rename nodes later according to class table!!!
dPSI <- dPSI_vt_full_chaging
dPSI$RANGE <- NULL
iis <- match(colnames(dPSI), class$Gene.Symbol)
iis

colnames(dPSI)[which(is.na(iis))]
# "AA1"       "AA2"       "AA3"       "AA4"       "AA5"       "AA6"       "AA7"      
# "AA8"       "AA9"       "C1orf55"   "C1orf55_b" "CWC22_b"   "HFM1_b"    "IKcon"    
#"PPIL1"     "PRPF8con"  "SF3B1con"  "SMU1con"   "SRPK2_b"   "SRRM4"     "SRRT"     
# "XAB2_b"
##########################################
#sdPSI
getwd()
write.table(round(dPSI, digits = 2), "ES_dPSI.tab",  sep = "\t")
##########################################
#single scaled
deltascaled <- scale(dPSI)# scaled by columns (KDs)
write.table(deltascaled, "ES_vt_sscaled.tab",  sep = "\t")
####################################################
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
write.table(dsscaled, "ES_vt_dscaled.tab",  sep = "\t")
####################################################


