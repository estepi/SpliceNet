library(pheatmap)
library(tidyr)
library(limma)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape) 
library(igraph)
##################################################################
setwd("~/Documents/SpliceNetRes/cor03/")
##############################
sampleClass<-read.csv("~/Documents/summaryLinks/class_colors.tab", sep="\t", header = T)
edgelists<-read.table("edgelist.tab", header = F, stringsAsFactors = F ,sep="\t")
edgelists[1:4,]; 
dim(edgelists)#2
finalDF<-data.frame()
linksList<-list()
length(linksList)
nodesList<-list()
for (i in 1:nrow(edgelists))
{
  file<-edgelists$V1[i]
  print(file)
  name<-gsub("_edgelist.tab", "", edgelists$V1[i])
  print(name)
  edgeListLC<-read.table(file, header = T, row.names = 1);
  print(head(edgeListLC));
  print(dim(edgeListLC))
  ###########################################################
    #Links
  linksList[[i]]<-edgeListLC$link
  names(linksList)[[i]]<-name
  ###########################################################
  #nodes
  n1<-unique(c(as.character(edgeListLC$source),
               as.character(edgeListLC$target)))
  nodesList[[i]]<-n1
  names(nodesList)[[i]]<-name
  }
##################################################################
head(nodesList)
length(nodesList)
length(linksList)

allNodes<-unique(unlist(nodesList)); length(allNodes)
allNodesDF<-data.frame(matrix(0, ncol=length(nodesList), 
                                  nrow=length(allNodes)), 
                                  row.names = allNodes)
#############################################
#select KD:
for (i in 1:length(nodesList))
{
  ii<-match(nodesList[[i]], rownames(allNodesDF))
  colnames(allNodesDF)[i]<-names(nodesList)[i]
  allNodesDF[ii,i]<-1  
}
colnames(allNodesDF)
head(allNodesDF)
##################################################################
#finalMatrixCOmparison
dim(allNodesDF)# 209
KDsSum<-rowSums(allNodesDF)
KDsSum
summary(KDsSum)
KDsSumByNetworks<-colSums(allNodesDF)
summary(KDsSumByNetworks)
sort(KDsSum)
#network de 100 nodos:
topNodes<-allNodesDF
dim(topNodes)
head(topNodes)
#esto se puede plotear como heatmap
ii<-match(rownames(topNodes), sampleClass$Gene.Symbol)

topNodes$order<-sampleClass$ORDER[ii]
topNodes<-topNodes[!is.na(ii),]
dim(topNodes)
head(topNodes)
topNodesOnto<-topNodes[order(topNodes$order),]
topNodesOnto[1:5,1:5]
topNodesOnto$order<-NULL
head(sampleClass)
pdf("presence_ausence_KDS.pdf",   width = 2,  height = 12)
pheatmap(topNodesOnto,
         fontsize_row = 2, 
         fontsize_col = 10,
        color =c("red","white"),
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        legend_breaks = 2)
dev.off()
#######################################
allLinks<-unique(unlist(linksList)); length(allLinks)  

#3159
allLinksDF<-data.frame(matrix(0, 
                              ncol=length(linksList), 
                              nrow=length(allLinks)), 
                              row.names = allLinks)
head(allLinksDF)
dim(allLinksDF)#2794
rownames(allLinksDF)[1:10]
head(allLinksDF)
rownames(allLinksDF)
#all possible conections:
for (i in 1:length(linksList))
  {
    print(i)  
    ii<-match(linksList[[i]], rownames(allLinksDF))
    
    colnames(allLinksDF)[i]<- names(linksList)[i]
    ff<-ii
    ff<-ff[!is.na(ff)]
    print(paste(ff,i))
    allLinksDF[ff,i]<-1  
    }
dim(allLinksDF)#1720
allLinksDF
head(allLinksDF)
###################################
colnames(allLinksDF)
###################################
totalLinksByNetwork<-colSums(allLinksDF)
summary(totalLinksByNetwork)
###################################
allLinksDF$occur<-rowSums(allLinksDF);head(allLinksDF)
allLinksDF<-allLinksDF[order(allLinksDF$occur, decreasing = T),]

allLinksDFTOP<-allLinksDF[order(allLinksDF$occur, decreasing = T)[1:500],]
allLinksDFLeast<-allLinksDF[order(allLinksDF$occur)[1:20],]
###################################
head(allLinksDF)
# 1 all no filter
linksPlot<-allLinksDFTOP
head(linksPlot)
########################################################
linksPlot$occur<-NULL
head(linksPlot)
colSums(linksPlot)
########################################################
pdf("score_linksTOP500.pdf",  width = 2,  height = 12)
    pheatmap(linksPlot,
             fontsize_row = 2, 
             fontsize_col = 10,
             color = colorRampPalette(c("white","red"))(2),
             cluster_cols = FALSE,
             cluster_rows = FALSE 
             )
    dev.off()
########################################################
ex1<-grep("SNRNP200", rownames(allLinksDF))
linksPlot<-allLinksDF[unique(ex1),]    
linksPlot
#######################################
linksPlot$occur<-NULL
head(linksPlot)
########################################################
pdf("score_links_BRR2.pdf",  width = 4,  height = 12)
  pheatmap(linksPlot,
           fontsize_row = 14, 
           fontsize_col = 14,
           color = colorRampPalette(c("white","red"))(2),
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           angle_col = "90"
  )
dev.off()

  ex1<-c(
    grep("IK", rownames(allLinksDF)),
    grep("CHERP", rownames(allLinksDF)),
    grep("SR140", rownames(allLinksDF)),
    grep("RBM17", rownames(allLinksDF)))
dim(allLinksDF)

ex1<-  grep("IK", rownames(allLinksDF))
linksPlot<-allLinksDF[unique(ex1),]    
linksPlot
linksPlot$occur<-NULL
head(linksPlot)
########################################################
pdf("score_links-IK.pdf",  width = 4,  height = 12)
pheatmap(linksPlot,
         fontsize_row = 14, 
         fontsize_col = 14,
         color = colorRampPalette(c("white","red"))(2),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         angle_col = "90"
)
dev.off()
