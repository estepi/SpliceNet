library(ggplot2)
library(igraph)
library(scales)
library(pheatmap)
library(plyr)
library(reshape) 
##################################################################
setwd("~/Documents/SpliceNetRes/GC/networks/cor03/")
rankings<-read.table("centrality.txt", header = F, stringsAsFactors = F ,sep="\t")
rankings
sampleClass<-read.csv("~/Documents/summaryLinks/class_colors.tab", sep="\t", header = T)
head(sampleClass)
##################################################################
                  metrics<-read.table(rankings$V1[1]);
                  head(metrics)
                  finalDF<-data.frame()
                  nodesList<-list()
                  nodesListNames<-list()
                  
                  for (i in 1:nrow(rankings))
                  {
                    file<-rankings$V1[i]
                    print(file)
                    name<-gsub("_CentralityRanking.tab", "", rankings$V1[i])
                    print(name)
                    metrics<-read.table(file);
                    ###########################################################
                    #nodes
                    head(metrics)
                    mdg<-metrics$PR
                    names(mdg)<-rownames(metrics)
                    nodesList[[i]]<-mdg
                    nodesListNames[[i]]<-rownames(metrics)
                    names(nodesList)[[i]]<-name
                  }
                  head(finalDF)
                  #total Degree: suma degree / 2
                  allNodes<-unique(unlist(nodesListNames)); length(allNodes)
                  allNodesDF<-data.frame(matrix(0, ncol=length(nodesList), nrow=length(allNodes)), row.names = allNodes)
                  head(allNodesDF)
                  dim(allNodesDF)#250 nodos
                  
                  for (i in 1:length(nodesList))
                  {
                    ii<-match(names(nodesList[[i]]), rownames(allNodesDF))
                    colnames(allNodesDF)[i]<-names(nodesList)[i]
                    allNodesDF[ii,i]<-as.numeric(nodesList[[i]]*100)
                  }
         
#####################################
head(allNodesDF)
forPlot<-allNodesDF
                  #paso a 0s y 1s
                  rownames(allNodesDF)                  
                  forPlot$gene<-rownames(forPlot)
                  forPlotMelt <- melt(forPlot, variable_name = "gene")
                  colnames(forPlotMelt)[2]<-"network"
                  ii<-match(forPlotMelt$gene, sampleClass$Gene.Symbol)
                  forPlotMelt$order<-sampleClass$ORDER[ii]
                  head(forPlotMelt)
                  forPlotMelt
                  head(forPlotMelt)
                  ##################################################################
                  forPlotMelt$class<-sampleClass$Class...family.small[ii]
                  forPlotMelt$color<-sampleClass$color[ii]
                  head(forPlotMelt)
                  table(forPlotMelt$class)
                  colnames(forPlotMelt)
                  #################################################
                  complex<-unique(forPlotMelt$class)
                  head(complex)
                  cc<-match(complex, sampleClass$Class...family.small)
                  color<-sampleClass$color[cc]
                  names(color)<-complex
                  color
                  head(forPlotMelt)
                  as.factor(forPlotMelt$network)
                  write.table(forPlotMelt, "summary_ggplot_PR.tab", sep="\t", col.names = NA)
                
filter<-forPlotMelt[forPlotMelt$order<50,]            
filter<-forPlotMelt
#################################################
pdf("SummaryPR_noLegend_all.pdf",
    width = 2,
    height = 12)
ggplot(data = filter) +
  geom_tile(aes(
    x = network,
    y = reorder(gene, order),
    fill = class,
    alpha = value
  )) +
  theme_classic() +
  scale_fill_manual(values = color) + #assign tissues colors
  scale_alpha_identity() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 10
  )) +
  theme(axis.text.y = element_text(size = 5)) +
  theme(legend.position = "none")
dev.off()

#################################################

