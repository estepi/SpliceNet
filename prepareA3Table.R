library(MatrixGenerics)
###########################################
dPSI_full = read.delim("~/Dropbox (CRG)/LAB_VALCARCEL/Network/tables/dPSI_full_No_Nas.txt",
                       dec = ",")
dim(dPSI_full)

table(dPSI_full$COMPLEX)
dPSI_full_alt3 = dPSI_full[dPSI_full$COMPLEX %in% "Alt3", ]
dim(dPSI_full_alt3)

dPSI_full_alt3[grep("1/1", rownames(dPSI_full_alt3)),]
mEv<-matrix(unlist(strsplit(rownames(dPSI_full_alt3), "-")), ncol=2, byrow = T)
head(mEv)
tt<-table(mEv[,1])
tt[tt>1]
hist(tt)
names(tt)[1]

mEv[duplicated(mEv[,1]),]
mEv["HsaALTA1005410",]

table(matrix(unlist(strsplit(rownames(dPSI_full_alt3), "-")), ncol=2, byrow = T)[,1])

summary(dPSI_full_alt3$LENGTH)
head(dPSI_full_alt3[,1:10])

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

mEv<-matrix(unlist(strsplit(rownames(dPSI_full_chaging), "-")), ncol=2, byrow = T)
head(mEv)
tt<-table(mEv[,1])
tt[tt>1]
hist(tt)
names(tt)[1]
ii<-grep("HsaALTA1038016", rownames(dPSI_full_chaging))

dPSI_full_chaging[ii, 1:10]

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
############################################
colnames(dPSI)
summary(dPSI$LENGTH)

############################################s
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
write.table(dsscaled, "A3_all_dscaled.tab",  sep = "\t")
write.table(deltascaled, "A3_all_sscaled.tab",  sep = "\t")
####################################################

####################################################

