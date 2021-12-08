library(ggplot2)
library(reshape) 
library(scales)
###############################################
output<-"~/Documents/SpliceNetRes/KDsExon10Res/"
#define where you want to put output
dir(output)
folder<-"~/Documents/SpliceNetRes/KDsExon10Res/"
setwd(folder)
featuresClass<-read.table("~/Documents/MOFA/clean/r3/matt/res/features_class_order.csv", sep="\t", header = T)
head(featuresClass)

#list of folders which
dir(folder)
exp<-data.frame( dir(folder,pattern="^exons_do_c"), stringsAsFactors = F)
colnames(exp)<-"exp"

featuresTable<-data.frame()
directionTable<-data.frame()
meanTable<-data.frame()
meanTableBg<-data.frame()
dir()
###############################################
for (i in 1:nrow(exp))
{
  dset<-exp$exp[i]
  name<-gsub("exons_do_c_", "", dset)
  print(name)
  #################
  paste(dset,"overview_of_features_and_comparisons.tab", sep="/")
  
  col<-"PVAL_MANN_WHITNEY_U_TEST_down_VS_bgConst"
  featuresTable<-rbind(featuresTable,
                       data.frame(t(read.table(file= 
                                              paste(dset,"overview_of_features_and_comparisons.tab", sep="/"),
                                               header = T, row.names = 1, sep="\t"))) [col,]  )
  ########################################################
  colDir<-"QUAL_MANN_WHITNEY_U_TEST_down_VS_bgConst"

  directionTable<-rbind(directionTable,
                        data.frame(t(read.table(file= 
                        paste(dset,"overview_of_features_and_comparisons.tab", sep="/"),
                        header = T, row.names = 1, sep="\t"))) [colDir,]  )
  rownames(featuresTable)[i]<-name
  rownames(directionTable)[i]<-name
  ###########################################################
  meanDirBg<-paste("bgConst", "MEAN", sep="_")
  print(meanDirBg)
  meanTableBg<-rbind(meanTableBg,
  data.frame(t(read.table(file= 
  paste(dset,"overview_of_features_and_comparisons.tab", sep="/"),
  header = T, row.names = 1, sep="\t"))) [meanDirBg,]  )
  
  rownames(meanTableBg)[i]<-name
  ##################################################
  meanDir<-"down_MEAN"
  print(meanDir)
  meanTable<-rbind(meanTable,
                   data.frame(t(read.table(file= 
                                             paste(dset,"overview_of_features_and_comparisons.tab", sep="/"),
                                           header = T, row.names = 1, sep="\t"))) [meanDir,]  )
  
  rownames(meanTable)[i]<-name
}
###########################################################
head(featuresTable); dim(featuresTable)
head(directionTable); dim(directionTable)
head(meanTable); dim(meanTable)
head(meanTableBg); dim(meanTableBg)

write.table(featuresTable, file="features_all.tab", sep="\t", col.names = NA)
write.table(directionTable, file="direction_all.tab", sep="\t", col.names = NA)
write.table(meanTable, file="meanTable_all.tab", sep="\t", col.names = NA)
write.table(meanTableBg, file="meanTable_BgC.tab", sep="\t", col.names = NA)
###########################################################
featuresTable<-read.table("features_all.tab", 
                          sep="\t", header = T, row.names = 1)
directionTable<-read.table("direction_all.tab",
                           sep="\t", header = T, row.names = 1)
meanTable<-read.table("mean_all.tab", 
                      sep="\t", header = T, row.names = 1)
###########################################################
Numeric<-data.frame(matrix(unlist(apply(meanTable,1,
                                        as.numeric)), 
                           ncol=ncol(meanTable),
                           nrow=nrow(meanTable), 
                           byrow = T))
head(Numeric)
colnames(Numeric)<-colnames(meanTable)
rownames(Numeric)<-rownames(meanTable)
dim(Numeric)
NumericSc<-data.frame(apply(Numeric,2,  function(x)
        { scales::rescale(x, to=c(0.2,1) )}  ))

colnames(NumericSc)<-colnames(meanTable)
rownames(NumericSc)<-rownames(meanTable)
dim(NumericSc)
###############################################
#featuresClass<-positionClass
head(featuresClass)
featuresClass<-featuresClass[order(featuresClass$ORDER),]
###########################################################
ii<-match(featuresClass$class, colnames(featuresTable))
figure<-featuresTable[,ii]

forPlot<-figure
forPlot$dset<-rownames(forPlot)
forPlotMelt <- melt(forPlot, 
                    id = "dset")
colnames(forPlotMelt)[2]<-"feature"
table(forPlotMelt$feature)
table(forPlotMelt$dset)

cc<-match(forPlotMelt$feature, featuresClass$class)

forPlotMelt$order<-featuresClass$ORDER[cc]
forPlotMelt$shortName<-featuresClass$FIGURE[cc]
head(forPlotMelt)
#############################################
#add direction:
head(directionTable)
dd<-match(featuresClass$class, colnames(directionTable))
figuredd<-directionTable[,dd]
head(figuredd)
class(figuredd)
forPlotdd<-figuredd
forPlotdd$dset<-rownames(forPlotdd)
forPlotMeltdd <- melt(forPlotdd, 
                    id = "dset")
colnames(forPlotMeltdd)[2]<-"feature"
dim(forPlotMelt)
dim(forPlotMeltdd)
###############################################
mm<-match(featuresClass$class, 
          colnames(NumericSc))
figuremm<-NumericSc[,mm]
head(figuremm)
class(figuremm)
forPlotmm<-figuremm
forPlotmm$dset<-rownames(forPlotmm)
forPlotMeltmm <- melt(forPlotmm, 
                      id = "dset")
colnames(forPlotMeltmm)[2]<-"feature"
dim(forPlotMeltmm)
dim(forPlotMeltdd)
###############################################
figureNum<-Numeric[,mm]
head(figuremm)
write.table(figureNum, "meanValues.tab", sep="\t")
forPlotVal<-figureNum
forPlotVal$dset<-rownames(forPlotVal)
forPlotMeltVal <- melt(forPlotVal, 
                      id = "dset")
colnames(forPlotMeltVal)[2]<-"feature"
dim(forPlotMeltVal)
dim(forPlotMeltdd)
head(forPlotMeltdd)
dim(forPlotMelt)

combined<-cbind(forPlotMelt, 
                forPlotMeltdd,
                forPlotMeltmm,
                forPlotMeltVal)
table(combined$feature)
table(combined$dset)

head(combined)
dim(combined)
combined$dset<-NULL
combined$feature<-NULL
colnames(combined)[6]<-"direction"
combined$direction[grep(">",combined$direction)]<-"l"
combined$direction[grep("<",combined$direction)]<-"g"
combined$direction[grep("=",combined$direction)]<-"eq"
table(combined$direction)
head(combined)
combined$dset<-NULL
combined$feature<-NULL
colnames(combined)[7]<-"meanVal"
combined$meanVal
summary(combined$meanVal)
combined$dset<-NULL
combined$feature<-NULL
colnames(combined)[1]<-"PVAL"
colnames(combined)[8]<-"MEAN"
###################################################
combined$meanVal[combined$direction=="eq"]<-0
combined$dset <- factor(combined$dset,
              levels=rev(unique(combined$dset)))
head(combined)
###################################
combined$MEAN<-round(combined$MEAN, digits = 2)
#########################################################
head(combined)  
combinedBg<-combined
#melt
meanBg<-as.numeric(meanTableBg[1,])
names(meanBg)<-colnames(meanTableBg)
meanBgF<-data.frame(matrix(ncol=ncol(combined), 
                    nrow=length(unique(combined$feature))))
colnames(meanBgF)<-colnames(combined)
meanBgF$dset<-"bg"
meanBgF$direction<-"eq"
meanBgF$PVAL<-1
meanBgF$feature<-unique(combined$feature)
meanBgF$order<-combined$order[match(meanBgF$feature, combined$feature)]
meanBgF$shortName<-combined$shortName[match(meanBgF$feature, combined$feature)]
meanBgF$meanVal<-0
meanBgF$MEAN<-round(meanBg[match(meanBgF$feature, names(meanBg))], digits = 2)
head(meanBgF)
head(combined)
combinedBG<-rbind(combined, meanBgF)
###############################################
combinedBG$dsetOrder<-rep(1, nrow(combinedBG))
levels(combinedBg$dset)
pdf("Exons_all_bg0_kdsRe.pdf",
    width=18, height=18)
ggplot(combinedBG, aes(y= reorder(dset,dsetOrder),
                       x=reorder(shortName,order))) + 
  geom_tile( aes( fill=direction, 
                  alpha=meanVal),
             color = "black", size=0)   + 
  geom_text(aes(label=MEAN), size=4, fontface='bold') +
  scale_alpha_identity() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0, 
                                   size=20))+
  theme(axis.text.y = element_text( size=10))+
  scale_fill_manual(values=c("white","red","blue")) +
  scale_x_discrete(position = 'top')+
  scale_y_discrete(position = 'right')
dev.off()
###############ordenar por la mean de 1 variable:
