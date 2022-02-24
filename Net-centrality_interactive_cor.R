##############################################################
library(parallel)##multicores
library(Hmisc)# correlation matrix
library(tibble)#rownmaes_to_columns
library(tidyr)#gather
library(dplyr)#left join
library(igraph)# networks
library(scales)#scaling
library(reshape) #for plot
library(ggplot2)#for plot
library(pheatmap)#for plot
#

t_Zvals<-read.table("/home/estepi/Documents/SpliceNetData/input/ES_subset3500_sscaled.tab")
sampleData <- as.matrix(t_Zvals)
sampleClass<-read.csv("~/Documents/summaryLinks/class_colors.tab", sep="\t", header = T)
name<-"ES"
start<- 0.1
end<- 0.6
interval<-0.02
ncores <- 1
fdr<-0.01
scripts<-"/home/estepi/Documents/SpliceNet"
sc3 <- paste(scripts, "functions_centrality_MC_cor.R", sep = "/")
#esta funcion se carga de sc3
source(sc3)
DGList <- getCentralityByCOR(sampleData, start, end, interval, fdr)
#########################################################
print("Finish computation. Lets plot")
DGdf <- DGList[[1]]
DGdf[is.na(DGdf)] <- 0

FileName<-paste(name, start, end, sep="-")
degreeFile <- paste(FileName, "degree.tab", sep = "_")
write.table(DGdf, degreeFile, sep = "\t", col.names = NA)
#########################################################
NormDegreeFile <- paste(FileName, "norm_degree.tab", sep = "_")
NDGdf <- DGList[[2]]
NDGdf[is.na(NDGdf)] <- 0
write.table(NDGdf, NormDegreeFile, sep = "\t", col.names = NA)
#########################################################
forPlot <- NDGdf
forPlot$gene <- rownames(forPlot)
forPlotMelt <- melt(forPlot, variable_name = "gene")
colnames(forPlotMelt)[2] <- "threshold"

ii <- match(forPlot$gene, sampleClass$gene.name.VT)
forPlotMelt$order <- sampleClass$ORDER[ii]
forPlotMelt$rescale <- rescale(forPlotMelt$value, to = c(0.1, 1))
ii <- match(forPlotMelt$gene, sampleClass$gene.name.VT)
forPlotMelt$class <- sampleClass$Class...family.small[ii]
forPlotMelt$color <- sampleClass$color[ii]
complex <- unique(forPlotMelt$class)
cc <- match(complex, sampleClass$Class...family.small)
color <- sampleClass$color[cc]
names(color) <- complex
#################################################
tilePlot <- paste(FileName, "NormDegree.pdf", sep = "_")
pdf(tilePlot, width = 4, height = 12)
ggplot(data = forPlotMelt) +
  geom_tile(aes(
    x = threshold,
    y = reorder(gene, order),
    fill = class,
    alpha = rescale
  )) +
  theme_classic() +
  scale_fill_manual(values = color) + #assign tissues colors
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 10
  )) +
  theme(axis.text.y = element_text(size = 3)) +
  theme(legend.position = "none")
dev.off()
#################################################
curvesPlot <- paste(FileName, "NormDegree_curves_order.pdf", sep = "_")
pdf(curvesPlot, width = 11.69,  height = 8.27)
ggplot(data = forPlotMelt,
       aes(
         x = factor(threshold , level = seq(start, end, interval)),
         y = value,
         group = order,
         colour = order
       )) +
  theme_minimal() +
  geom_line()
dev.off()
#################################################
curvesPlotClass <- paste(FileName, "NormDegree_curves_class.pdf", sep = "_")
pdf(curvesPlotClass, width = 11.69,  height = 8.27)
ggplot(data = forPlotMelt,
       aes(
         x = factor(threshold , level = seq(start, end, interval)),
         y = value,
         group = order,
         colour = class
       )) +
  theme_minimal() +
  geom_line()
dev.off()
##############################################
#cluter factors (aacording norm degree)
#add annotation
#correct scale
###############################################
paletteLength=100
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
ClassID<-data.frame(gene=sampleClass$gene.name.VT,
                    class=sampleClass$Class...family.small)
rownames(ClassID)<-sampleClass$gene.name.VT
ClassID$gene<-NULL
ClassColors<-data.frame(class=unique(sampleClass$Class...family.small))
ClassColors$RGBcolor <- sampleClass$color[match(ClassColors$class,
                                                sampleClass$Class...family.small)]

ClassColorsVector<-ClassColors$RGBcolor
names(ClassColorsVector)<-ClassColors$class
ann_colors = list(class =  ClassColorsVector)

heatmapPDF<-paste(FileName, "ND_heatmap.pdf",sep="_")
pdf(heatmapPDF, width =6,  height = 20)
pheatmap(
  NDGdf,
  fontsize_col = 10,
  fontsize_row = 6,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = myColor,
  annotation_row =   ClassID,
  annotation_colors = ann_colors)
dev.off()

heatmapPDF<-paste(FileName, "DG_heatmap.pdf",sep="_")
pdf(heatmapPDF, width =6,  height = 20)
pheatmap(
  DGdf,
  fontsize_col = 10,
  fontsize_row = 6,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = myColor,
  annotation_row =   ClassID,
  annotation_colors = ann_colors)
dev.off()