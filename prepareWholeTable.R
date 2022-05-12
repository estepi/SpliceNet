library(MatrixGenerics)
###########################################
dPSI_full = read.delim("~/Documents/SpliceNetData/data/dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)

MIN <- rowMins(as.matrix(dPSI_full[, 18:322]))
MAX <- rowMaxs(as.matrix(dPSI_full[, 18:322]))
RANGE = MAX - MIN
Sds <- rowSds(as.matrix(dPSI_full[, 18:322]))
dPSI_full = cbind(dPSI_full, RANGE = RANGE)
dPSI_full = cbind(dPSI_full, Sds = Sds)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_full_chaging = dPSI_full[dPSI_full$RANGE > 5, ]
dim(dPSI_full_chaging) #41870   324
table(dPSI_full_chaging$COMPLEX)

#Alt3  Alt5   ANN    C1    C2    C3    IR   MIC     S 
#6093  3676  1088  1513  1359  1772 15380    72 10917 

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
colnames(dPSI_full_chaging)
#original table is names according to Gene.Symbol
class <-
  read.delim(
    "~/Documents/SpliceNetData/data//class_colors_2020.txt",
    header = T
  )
#so for the network I rename nodes later according to class table!!!
colnames(dPSI_full_chaging)[!colnames(dPSI_full_chaging) %in% class$Gene.Symbol]
dPSI <- dPSI_full_chaging
dPSI <- dPSI[, 18:ncol(dPSI)]
dPSI$RANGE <- NULL
dPSI$Sds <- NULL
colnames(dPSI)[!colnames(dPSI) %in% class$Gene.Symbol]
dim(dPSI)#305
iis <- match(colnames(dPSI), class$Gene.Symbol)
iis[is.na(iis)]
colnames(dPSI) <- class$gene.name.VT[iis]
##########################################
dPSI_print<-cbind(dPSI,
                  AVcontrols=dPSI_full_chaging$AVcontrols)


deltascaled <- scale(dPSI)# scaled by columns (KDs)

dsscaled <- t(scale(t(deltascaled)))    #scaled by events
#double scaled
#single scaled
####################################################
setwd("~/Documents/SpliceNetData/data/")
dim(dPSI_print)

write.table(dPSI_print, "dPSI_all.tab",  sep = "\t")
write.table(dsscaled, "all_dscaled.tab",  sep = "\t")
write.table(deltascaled, "all_sscaled.tab",  sep = "\t")
####################################################


