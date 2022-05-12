library(MatrixGenerics)
###########################################
dPSI_full = read.delim("../SpliceNetData//dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)
table(dPSI_full$COMPLEX)
dPSI_full_alt3 = dPSI_full[dPSI_full$COMPLEX %in% "Alt3", ]

dPSI_full_alt3[grep("1/1", rownames(dPSI_full_alt3)),]
mEv<-matrix(unlist(strsplit(rownames(dPSI_full_alt3), "-")), ncol=2, byrow = T)

MIN <- rowMins(as.matrix(dPSI_full_alt3[, 18:322]))
MAX <- rowMaxs(as.matrix(dPSI_full_alt3[, 18:322]))
RANGE = MAX - MIN
Sds <- rowSds(as.matrix(dPSI_full_alt3[, 18:322]))
dPSI_full_alt3 = cbind(dPSI_full_alt3, RANGE = RANGE)
dPSI_full_alt3 = cbind(dPSI_full_alt3, Sds = Sds)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_full_chaging = dPSI_full_alt3[dPSI_full_alt3$RANGE > 5, ]
dim(dPSI_full_chaging) #6093 rows
table(dPSI_full_chaging$COMPLEX)

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
############################################
dPSI <- dPSI[, 18:ncol(dPSI)]
dPSI$RANGE <- NULL
dPSI$Sds <- NULL
iis <- match(colnames(dPSI), class$Gene.Symbol)
colnames(dPSI) <- class$gene.name.VT[iis]
##########################################
setwd("../SpliceNetData/")
##########################################
#sdPSI
write.table(round(dPSI, digits = 2), "A3_dPSI.tab",  sep = "\t")
##########################################
#single scaled
deltascaled <- scale(dPSI)# scaled by columns (KDs)
write.table(deltascaled, "A3_all_sscaled.tab",  sep = "\t")
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
deltascaled[1:5,300:305]
write.table(dsscaled, "A3_all_dscaled.tab",  sep = "\t")
####################################################
