library(MatrixGenerics)
###########################################
dPSI_full = read.delim("../SpliceNetData//dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)

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
write.table(round(dPSI, digits = 2), "A5_dPSI.tab",  sep = "\t")
##########################################
#single scaled
deltascaled <- scale(dPSI)# scaled by columns (KDs)
write.table(deltascaled, "A5_all_sscaled.tab",  sep = "\t")
####################################################
#double scaled
dsscaled <- t(scale(t(deltascaled)))    #scaled by events
write.table(dsscaled, "A5_all_dscaled.tab",  sep = "\t")
####################################################


