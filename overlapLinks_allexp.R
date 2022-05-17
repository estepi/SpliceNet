library(ggplot2)
library(UpSetR)
library(igraph)
library(ggpubr)
##################################################################
setwd("SpliceNetRes")
class<-read.delim("../SpliceNetData/class_colors_2020.txt", header = T)
##################################################################
expected <- read.table(
  "expected_alphSort.tab",
  sep = "\t",
  header = T,
  row.names = 1
)

checks <-
  which(is.na(match(expected$source, class$gene.name.VT)))
length(checks)
checkt <-
  which(is.na(match(expected$target, class$gene.name.VT)))
length(checkt)
remove <- union(checks, checkt)
length(remove)#786
#nodos que no estan en las redes pero estan en las DBs pues son interesantes
expected <- expected[-remove, ]
dim(expected)
length(unique(expected$l11))
class(expected)
#load libraries
##############################
sourcedf <- as.data.frame(table(expected$short))
colnames(sourcedf) <- c("source", "times")
head(sourcedf)
#figure 1
##################################################################
colnames(expected)
listSources <- list()
nodesSourceList <- list()
expectedTop <- expected
sources <- names(table(expectedTop$short))
length(unique(expectedTop$l11))
length(unique(expected$l11))

for (i in 1:length(sources))
{
  print(sources)[i]
  listSources[[i]]  <- expected$l11[expected$short == sources[i]]
  names(listSources)[i] <- sources[i]
  #######################################
  n1 <-
    unique(c(
      as.character(expected$source[expected$short == sources[i]]),
      as.character(expected$target[expected$short == sources[i]])
    ))
  nodesSourceList[[i]] <- n1
  names(nodesSourceList)[[i]] <- sources[i]
}
##################################################################
#evidence non from RNA-Seq
NoRNASeq<-listSources
names(NoRNASeq)[4]
names(NoRNASeq)[1]
NoRNASeq[1]<-NULL
NoRNASeq[1]<-NULL
NoRNASeq[1]<-NULL
names(NoRNASeq)[1]
#repetir varias veces!!!
length(NoRNASeq)#8
names(NoRNASeq)
#exclude RNA-Seq, FIG 3
#se podria ordenar mejor los sets
RNASeq<-listSources
names(RNASeq)
names(RNASeq)[4]
RNASeq[4:6]<-NULL
RNASeq[5:8]<-NULL
length(RNASeq)#4
names(RNASeq)#OK
#######################################
#FIG 4: RNA-Seq cross vs LABCHIP
onlyRNASeq<-RNASeq
names(onlyRNASeq)
onlyRNASeq[4]<-NULL
onlyRNASeqLinks<-unlist(onlyRNASeq)
##################################################
#convierto en la lista de redes:
##################################################
edgelists <-
  read.table(
    "edgelist.tab",
    header = F,
    stringsAsFactors = F ,
    sep = "\t"
  )
edgelists
#agrego las top 500:
finalDF <- data.frame()
linksList <- list()
length(linksList)
nodesList <- list()
dim(edgelists)#10
edgelists

#All de caada network:
for (i in 1:nrow(edgelists)) {
  file <- edgelists$V1[i]
  print(file)
  name <- gsub("_edgelist.tab", "", edgelists$V1[i])
  print(name)
  edgeListLC <-
    read.table(edgelists$V1[i], header = T, row.names = 1)
  ###########################################################
  #Links hay que hacer el sort
  edgeListLCS <- data.frame(edgeListLC$source, edgeListLC$target)
  edgeListLCSort <-
    data.frame(t(unlist(apply(edgeListLCS[, 1:2], 1, function(x) {
      sort(x)
    }))))
  colnames(edgeListLCSort) <- c("source", "target")
  l11 <-
    paste(edgeListLCSort$source, edgeListLCSort$target, sep = "-")
  linksList[[i]] <- l11
  names(linksList)[[i]] <- name
  ###########################################################
  #nodes
  n1 <- unique(c(
    as.character(edgeListLC$source),
    as.character(edgeListLC$target)
  ))
  nodesList[[i]] <- n1
  names(nodesList)[[i]] <- name
}
##################################################################
lapply(linksList, length)
names(linksList)
linksList
###############################################
names(RNASeq[4])
linksListAll <- c(RNASeq[4], linksList)
names(linksListAll)
lapply(linksListAll, length)
getwd()
length(linksListAll)
lapply(linksListAll, length)
names(linksListAll)
upset(fromList(linksListAll),
      sets = c("Labchip",  "ES"))
#################################################################
#cual de las networks recupera mejor los datos experimentales?
names(NoRNASeq)
experimental <- NoRNASeq
experimental[[4]] <- NULL
allExp <- list(unique(unlist(experimental)))
length(allExp[[1]])

rnaseq_exp <- c(linksListAll, allExp)
length(rnaseq_exp)#6
names(rnaseq_exp)[6] <- "allExp"
names(rnaseq_exp)
lapply(rnaseq_exp, length)
#############################################################
#de acuerdo a esto, numeros:
upset(fromList(rnaseq_exp), sets = c("Labchip",  "allExp"))#237/541*100
upset(fromList(rnaseq_exp), sets = c("ES",  "allExp"))
upset(fromList(rnaseq_exp), sets = c("IR",  "allExp"))
upset(fromList(rnaseq_exp), sets = c("A5",  "allExp"))
upset(fromList(rnaseq_exp), sets = c("A3",  "allExp"))
#########################################
gL <-
  graph_from_data_frame(data.frame(matrix(
    unlist(strsplit(rnaseq_exp$Labchip, "-")), ncol = 2, byrow = T
  )))
gE <-
  graph_from_data_frame(data.frame(matrix(
    unlist(strsplit(rnaseq_exp$allExp, "-")), ncol = 2, byrow = T
  )))
intersection(gL, gE, keep.all.vertices = F)#237 links in 166 KDs
names(rnaseq_exp)

resultL <- data.frame(matrix(nrow = length(rnaseq_exp)) ,
                      row.names = names(rnaseq_exp))
resultE <- data.frame(matrix(nrow = length(rnaseq_exp)) ,
                      row.names = names(rnaseq_exp))


for(i in 1:length(rnaseq_exp)) {
  g2 <-
    graph_from_data_frame(data.frame(matrix(
      unlist(strsplit(rnaseq_exp[[i]], "-")), ncol = 2, byrow = T
    )))
  #########################################################################
  #against labchip
  resultL[i, 1] <- length(E(gL))
  resultL[i, 2] <- length(V(gL))
  resultL[i, 3] <- length(E(g2))
  resultL[i, 4] <- length(V(g2))
  resultL[i, 5] <-
    length(E(intersection(gL, g2, keep.all.vertices = F)))
  resultL[i, 6] <-
    length(V(intersection(gL, g2, keep.all.vertices = F)))
  resultL[i, 7] <- resultL[i, 5] / resultL[i, 3] * 100
  resultL[i, 8] <- resultL[i, 6] / resultL[i, 4] * 100
  #########################################################################
  #against PPI
  resultE[i, 1] <- length(E(gE))
  resultE[i, 2] <- length(V(gE))
  resultE[i, 3] <- length(E(g2))
  resultE[i, 4] <- length(V(g2))
  resultE[i, 5] <-
    length(E(intersection(gE, g2, keep.all.vertices = F)))
  resultE[i, 6] <-
    length(V(intersection(gE, g2, keep.all.vertices = F)))
  resultE[i, 7] <- resultE[i, 5] / resultE[i, 3] * 100
  resultE[i, 8] <- resultE[i, 6] / resultE[i, 4] * 100
  #########################################################################
}

colnames(resultL) <-
  c("benchL",
    "benchV",
    "g2L",
    "g2V",
    "intL",
    "intV",
    "intLp",
    "intVp")
head(resultL)
resultL$dset <- rownames(resultL)
write.table(resultL,
            file = "result_overlap_labchip_subset3500.tab",
            sep = "\t",
            col.names = NA)
##################################################
colnames(resultE) <-
  c("benchL",
    "benchV",
    "g2L",
    "g2V",
    "intL",
    "intV",
    "intLp",
    "intVp")
write.table(resultE,
            file = "result_overlap_Exp_subset3500.tab",
            sep = "\t",
            col.names = NA)
