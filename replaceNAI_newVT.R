replaceNAI<-function(INCFile,
                    OUTFile1,
                    OUTFile2
  ){
  #######################################################################
  #READ INPUT FILE
  INCData <- read.csv(INCFile, 
                      header=T, 
                      sep="\t", 
                      stringsAsFactors = F, row.names =2 )
  #######################################################################
  #subset PSI values
  INCData<-infile; 
  dim(INCData)
  rownames(infile)[1:10]
  colnames(infile)[1:6]
  table(infile$COMPLEX)
  onlyValues<-INCData[,7:ncol(INCData)]
  #subset QUAL values
  Qvals <- grep("\\.Q$", colnames(onlyValues))
  qualValues <- onlyValues[, Qvals]
  psiValues  <- onlyValues[,-Qvals]
  finalValues<-matrix(0, ncol=ncol(psiValues), nrow=nrow(psiValues))
  colnames(finalValues)<-colnames(onlyValues)[-Qvals]
  dim(finalValues)
  #######################################################################
#count Ns in QUAL values
Ns_comma<-t(apply(qualValues, 1, function(x){ str_count(x, "N,")}))
#######################################################################
#replace 
indexes_N<-which(Ns_comma !=0, arr.in=TRUE)
indexes_NewN3<-which(Ns_comma ==3, arr.in=TRUE);
finalValues[indexes_N]<-Ns_comma[indexes_N];head(finalValues)
correctedPsiVal<-psiValues;
correctedPsiVal[indexes_NewN3]<-"NAnew3"
originalNA<-is.na(psiValues)
correctedPsiVal[originalNA]<-"NAold"
dim(correctedPsiVal)
INCData[1:5,1:6]
fullNAS<-cbind(INCData[, 1:6], correctedPsiVal)
table(fullNAS$COMPLEX)
#####################################################
#compute and print STATS:
tVal<-length(originalNA)#
oNA<-length(which(originalNA))#
new3<-length(indexes_NewN3)/2# 
oNA2<-length(which(correctedPsiVal=="NAold"))#
new3b<-length(which(correctedPsiVal=="NAnew3"))#
m1<-paste(round(oNA/tVal*100),"%")
m2<-paste(round(new3/tVal*100),"%")
m3<-paste(round(new3b/tVal*100),"%")
message(paste("Num total PSI values:",tVal ) ) 
message(paste("NA original values:",  m1  ))
message(paste("NewN3 values:", m3 ))
message(paste("Total NA + NewN3 values:",m2) ) 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    #####################################################
#PRINT OUTPUT FILE
write.table(fullNAS,OUTFile1, sep="\t", col.names = NA )
message(paste("Final corrected table only NewNAs is written in:",OUTFile1) ) 
#####################################################
#correct IR:
  message("Number of events by type: ")
  print(table(fullNAS$COMPLEX) )   
  rownames(qualValues)[1:10]
  rownames(fullNAS)[1:10]
  fullNAS$COMPLEX[which(is.na(fullNAS$COMPLEX))]<-"S"
  ii<-as.character(fullNAS$COMPLEX)=="IR"
  table(fullNAS$COMPLEX)
  table(as.character(fullNAS$COMPLEX))
  length(which(is.na(ii)))
  sum(table(ii))
  length(ii)
  dim(qualValues)
  dim(fullNAS)
  colnames(qualValues)
  length(which(rownames(fullNAS)==rownames(qualValues)))
  head(fullNAS[fullNAS$COMPLEX=="IR",])
  IRCqual<-qualValues[fullNAS$COMPLEX=="IR",]
  head(IRCqual)
  rownames(IRCqual)<-rownames(qualValues)[fullNAS$COMPLEX=="IR"]
  IRCpsi<-fullNAS[fullNAS$COMPLEX=="IR" ,7:ncol(fullNAS)]
  head(IRCqual)
  #####################################################################
  pval<-apply(IRCqual, 2, function(x) 
    { 
      as.numeric(str_match(
        matrix( unlist(strsplit(x, ",")), ncol=6, byrow = T)[,5],"(.*)@")[,2])
    })
    warnings()
    pvalCorrected<-apply(pval, 2, p.adjust)
    pvalCorrected[pvalCorrected<0.05]<-"NAI"
    ######################################################################
    Nis<-t(apply(pvalCorrected, 1, function(x)
    { str_count(x, "NAI")}))
    indexes_NAI<-which(Nis !=0, arr.in=TRUE)
    #####################################################################
    dim(IRCpsi); dim(IRCqual)
    rownames(IRCpsi)[1:10]
    rownames(IRCqual)[1:10]
    correctedIRpsi<-IRCpsi; 
    rownames(IRCpsi)[1:10]
    #compute and print STATS IR
    tIR<-dim(IRCpsi)[1]*dim(IRCpsi)[2]
    tNewNAI<-dim(indexes_NAI)[1]
    
    m4<-paste(round(tNewNAI/tIR*100),"%")
  message(paste("Num total PSI values:",tVal ) ) 
  
  correctedIRpsi[indexes_NAI]<-"NAI"
  df1<-correctedPsiVal;   dim(df1)
  rownames(df1)[1:10]#OK
  
  df2<-correctedIRpsi;dim(df2)
  rownames(df2)[1:10]
  
  df1[match(rownames(df2) , rownames(df1)), ] <- df2
  fullNASIR<-cbind(INCData[, 1:6], df1)
  
  message("pval analyzed (for correction): ",length(pvalCorrected))
  message("pval corrected <0.05: ",m4)
  
  write.table(fullNASIR,OUTFile2, sep="\t", col.names = NA )
  message(paste("Final corrected table NewNAs+NAI is written in:",OUTFile2) ) 
  
  
}