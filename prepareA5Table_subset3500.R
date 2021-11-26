library(MatrixGenerics)
###########################################
dPSI_full = read.delim("~/Dropbox (CRG)/LAB_VALCARCEL/Network/tables/dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)

table(dPSI_full$COMPLEX)
dPSI_full_alt5 = dPSI_full[dPSI_full$COMPLEX %in% "Alt5", ]
dim(dPSI_full_alt5)
MIN <- rowMins(as.matrix(dPSI_full_alt5[, 18:322]))
MAX <- rowMaxs(as.matrix(dPSI_full_alt5[, 18:322]))
RANGE = MAX - MIN
Sds <- rowSds(as.matrix(dPSI_full_alt5[, 18:322]))
dPSI_full_alt5 = cbind(dPSI_full_alt5, RANGE = RANGE)
dPSI_full_alt5 = cbind(dPSI_full_alt5, Sds = Sds)
###### remove non changing rows I try different cut offs bigger then 0,5,10,15….
dPSI_full_chaging = dPSI_full_alt5[dPSI_full_alt5$RANGE > 5, ]
dim(dPSI_full_chaging) #3676 rows
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
#subset to 3500
subset<-dPSI[sample(1:nrow(dPSI), 3500,replace = FALSE),]
##########################################
#single scaled
deltascaled <- scale(subset)# scaled by columns (KDs)
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
####################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/standard/diffRho/subset3500/")
write.table(dsscaled, "A5_subset3500_dscaled.tab",  sep = "\t")
write.table(deltascaled, "A5_subset3500_sscaled.tab",  sep = "\t")
####################################################


