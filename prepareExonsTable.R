library(MatrixGenerics)
###########################################
dPSI_full = read.delim("~/Dropbox (CRG)/LAB_VALCARCEL/Network/tables/dPSI_full_No_Nas.txt",
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
colnames(dPSI_full_chaging)
#original table is names according to Gene.Symbol
class <-
  read.delim(
    "~/Dropbox (CRG ADV)/Personal_Estefania/Network/summaryLinks/edgelist/class_colors_2020.txt",
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
colnames(dPSI)

dim(dPSI)
colnames(dPSI)
##########################################
deltascaled <- scale(dPSI)# scaled by columns (KDs)
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
#double scaled
#single scaled
####################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/standard/diffRho/")
write.table(dsscaled, "ES_all_dscaled.tab",  sep = "\t")
write.table(deltascaled, "ES_all_sscaled.tab",  sep = "\t")
####################################################


