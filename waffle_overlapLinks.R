library(ggplot2)
library(igraph)
library(scales)
library(pheatmap)
library(plyr)
library(reshape) 
##################################################################
setwd("~/Documents/SpliceNetRes/all/cor03/")
sampleClass<-read.csv("~/Documents/summaryLinks/class_colors.tab", sep="\t", header = T)
head(sampleClass)
##############################
edgelists<-read.table("edgelist.tab", header = F, stringsAsFactors = F ,sep="\t")##################################################################
edgelist<-read.table(edgelists$V1[1]);
head(edgelist)
finalDF<-data.frame()
dim(edgelists)#2
linksList<-list()
linkListNames<-list()
                
for (i in 1:nrow(edgelists))
{
  file <- edgelists$V1[i]
  print(file)
  name <-    gsub("_edgelist.tab", "", edgelists$V1[i])
  print(name)
  edgelist <- read.table(file, header = T)
  #links
  head(edgelist)
  absCor <- edgelist$absCor
  names(absCor) <- edgelist$link
  linksList[[i]] <- absCor
  names(linksList)[[i]] <- name
  linkListNames[[i]] <- edgelist$link
  names(linkListNames)[[i]] <- name
}

linksList
head(finalDF)
#total Degree: suma degree / 2
allLinks <-   unique(unlist(linkListNames))
head(allLinks)
length(allLinks)#2794
allLinksDF <-
data.frame(matrix(0, 
                    ncol = length(linksList),
                    nrow = length(allLinks)), 
                    row.names = allLinks)
head(allLinksDF)
dim(allLinksDF)#2794

for (i in 1:length(linksList))
{
  names(linksList[[i]])
  rownames(allLinksDF)
  
  ii <- match(
  names(linksList[[i]]), rownames(allLinksDF)    )
  
  colnames(allLinksDF)[i] <- names(linksList)[i]
  allLinksDF[ii, i] <-  linksList[[i]]
  
}
#####################################
#remove NAs
forPlot<-allLinksDF
freq<-rowSums(allLinksDF>0)
table(freq)
top<-names(freq)[freq>=3]
length(top)#386
forPlot<-allLinksDF[top,]
dim(forPlot)
head(forPlot)
MeanCor<-rowMeans(forPlot)
head(MeanCor)
#paso a 0s y 1s
forPlot$link <- rownames(forPlot)
head(forPlot)

forPlotMelt <-
  melt(forPlot, variable_name = "link")
colnames(forPlotMelt)[2] <- "network"
##################################################################
mm<-match(forPlotMelt$link, names(MeanCor))
forPlotMelt$order<-MeanCor[mm]
filter<-forPlotMelt
head(filter)
#################################################
pdf("SummaryPR_noLegend_freq.pdf",
    width = 2,
    height = 12)
ggplot(data = filter) +
  geom_tile(aes(
    x = network,
    y = reorder(link, order),
    fill = "red",
    alpha = value
  )) +
  theme_classic() +
  scale_fill_manual(values = "red") + 
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    #hjust = 400,
    size = 10
  )) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 2)) +
  theme(legend.position = "none")
dev.off()

#####################################################
ik<-  grep("IK", filter$link)
linksPlot<-filter[unique(ik),]    

pdf("SummaryCor_noLegend_IK.pdf",
    width = 2,
    height = 2)
ggplot(data = linksPlot) +
  geom_tile(aes(
    x = network,
    y = reorder(link, order),
    fill = "red",
    alpha = value
  )) +
  theme_classic() +
  scale_fill_manual(values = "red") + 
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    #hjust = 400,
    size = 10
  )) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 5)) +
  theme(legend.position = "none")
dev.off()
####################################################
BRR2<-grep("SNRNP200", filter$link)
linksPlot<-filter[unique(BRR2),]    
dim(linksPlot)

pdf("SummaryCor_noLegend_BRR2.pdf",
    width = 2,
    height = 4)
ggplot(data = linksPlot) +
  geom_tile(aes(
    x = network,
    y = reorder(link, order),
    fill = "red",
    alpha = value
  )) +
  theme_classic() +
  scale_fill_manual(values = "red") + 
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    #hjust = 400,
    size = 10
  )) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 5)) +
  theme(legend.position = "none")
dev.off()

###########################################
CRNKL1<-grep("CRNKL1", filter$link)
linksPlot<-filter[unique(CRNKL1),]    
dim(linksPlot)
pdf("SummaryCor_noLegend_CRNKL1.pdf",
    width = 2,
    height = 4)
ggplot(data = linksPlot) +
  geom_tile(aes(
    x = network,
    y = reorder(link, order),
    fill = "red",
    alpha = value
  )) +
  theme_classic() +
  scale_fill_manual(values = "red") + 
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    #hjust = 400,
    size = 10
  )) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 5)) +
  theme(legend.position = "none")
dev.off()

