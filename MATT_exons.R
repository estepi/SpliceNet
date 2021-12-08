library(MatrixGenerics)
###########################################
dPSI_full = read.delim("~/Documents/SpliceNetData/dPSI_full_No_Nas.txt",
                                              dec = ",")
dim(dPSI_full)

dPSI_full_exons = dPSI_full[dPSI_full$COMPLEX %in% c("S", "C1", "C2", "C3", "ANN", "MIC"), ]
dim(dPSI_full_exons)

dPSI_full_exons[1:5,1:10]
getwd()
write.table(dPSI_full_exons, "../SpliceNetData/dPSI_full_exons.tab", sep="\t", quote = FALSE)
colnames(dPSI_full_exons)
summary(dPSI_full_exons$SDcontrols)

UP = dPSI_full_exons[dPSI_full_exons$AC008073.5 > 10, "EVENT"]
#altrernativ no changing
bg1 = dPSI_full_exons[abs(dPSI_full_exons$AC008073.5) < 5  &
                        dPSI_full_exons$AVcontrols > 5 &
                        dPSI_full_exons$AVcontrols < 95,  "EVENT"]
length(bg1)#500
#constituve not changing
bg0 = dPSI_full_exons[abs(dPSI_full_exons$AC008073.5) < 5 
                      &   dPSI_full_exons$AVcontrols < 5 |
                        dPSI_full_exons$AVcontrols > 95,  "EVENT"]
length(bg0)#139
