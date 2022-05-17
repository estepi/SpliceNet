library(ggplot2)
library(igraph)
library(scales)
library(pheatmap)
library(plyr)
library(reshape) 
##################################################################
setwd("SpliceNetData")
sampleClass <- read.csv("class_colors2020.tab", sep = "\t", header = T)
head(sampleClass)
##############################
edgelists <-
  read.table(
    "edgelists.txt",
    header = F,
    stringsAsFactors = F ,
    sep = "\t"
  )
##################################################################
edgelist <- read.table(edgelists$V1[1])

head(edgelist)
finalDF <- data.frame()
dim(edgelists)#2
linksList <- list()
linkListNames <- list()

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

#total Degree: suma degree / 2
allLinks <-   unique(unlist(linkListNames))
allLinksDF <-
  data.frame(matrix(0,
                    ncol = length(linksList),
                    nrow = length(allLinks)),
             row.names = allLinks)

for (i in 1:length(linksList))
{
  names(linksList[[i]])
  rownames(allLinksDF)
  
  ii <- match(names(linksList[[i]]), rownames(allLinksDF))
  
  colnames(allLinksDF)[i] <- names(linksList)[i]
  allLinksDF[ii, i] <-  linksList[[i]]
  
}
#####################################
#remove NAs
forPlot <- allLinksDF
freq <- rowSums(allLinksDF > 0)
top <- names(freq)[freq == 2]
forPlot <- allLinksDF[top, ]
MeanCor <- rowMeans(forPlot)
#paso a 0s y 1s
forPlot$link <- rownames(forPlot)

forPlotMelt <-
  melt(forPlot, variable_name = "link")
colnames(forPlotMelt)[2] <- "network"
##################################################################
mm <- match(forPlotMelt$link, names(MeanCor))
forPlotMelt$order <- MeanCor[mm]
filter <- forPlotMelt
#################################################
pdf("SummaryPR_noLegend_freq.pdf",
    width = 1,
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
  theme(axis.text.x = element_text(angle = 90,
                                   #hjust = 400,
                                   size = 10)) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(size = 2)) +
  theme(legend.position = "none")
dev.off()
#####################################################
