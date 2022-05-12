library(MatrixGenerics)
###########################################
dPSI_full = read.delim("../SpliceNetData//dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)

dPSI_full_exons = dPSI_full[dPSI_full$COMPLEX %in% c("S", "C1", "C2", "C3", "ANN", "MIC"), ]
dim(dPSI_full_exons)
MIN <- rowMins(as.matrix(dPSI_full_exons[, 18:322]))
MAX <- rowMaxs(as.matrix(dPSI_full_exons[, 18:322]))
RANGE = MAX - MIN
Sds <- rowSds(as.matrix(dPSI_full_exons[, 18:322]))
dPSI_full_exons = cbind(dPSI_full_exons, RANGE = RANGE)
dPSI_full_exons = cbind(dPSI_full_exons, Sds = Sds)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_full_chaging = dPSI_full_exons[dPSI_full_exons$RANGE > 5, ]
dim(dPSI_full_chaging) #16721 rows
table(dPSI_full_chaging$COMPLEX)
#ANN    C1    C2    C3   MIC     S 
#1088  1513  1359  1772    72 10917 

colnames(dPSI_full_chaging)[colnames(dPSI_full_chaging) == "LENG1_b"] <-
  "LENG1"
colnames(dPSI_full_chaging)[colnames(dPSI_full_chaging) == "RBM17con"] <-
  "RBM17"
colnames(dPSI_full_chaging)[colnames(dPSI_full_chaging) == "HFM1_b"] <-
  "HFM1"
colnames(dPSI_full_chaging)[colnames(dPSI_full_chaging) == "CCDC12_b"] <-
  "CCDC12"
colnames(dPSI_full_chaging)[colnames(dPSI_full_chaging) == "CDC5L_b"] <-
  "CDC5L"
dim(dPSI_full_chaging)
#original table is names according to Gene.Symbol
class <-
  read.delim(
    "../SpliceNetData//class_colors_2020.txt",
    header = T
  )
#so for the network I rename nodes later according to class table!!!
dPSI <- dPSI_full_chaging
dPSI <- dPSI[, 18:ncol(dPSI)]
dPSI$RANGE <- NULL
dPSI$Sds <- NULL
iis <- match(colnames(dPSI), class$Gene.Symbol)
colnames(dPSI) <- class$gene.name.VT[iis]
##########################################
setwd("../SpliceNetData/")
##########################################
#sdPSI
write.table(round(dPSI, digits = 2), "ES_dPSI.tab",  sep = "\t")
##########################################
#single scaled
deltascaled <- scale(dPSI)# scaled by columns (KDs)
write.table(deltascaled, "ES_all_sscaled.tab",  sep = "\t")
####################################################
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
write.table(dsscaled, "ES_all_dscaled.tab",  sep = "\t")
####################################################


