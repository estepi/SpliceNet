library(ggplot2)
library(UpSetR)
library(igraph)
library(ggpubr)
##################################################################
setwd("~/Dropbox (CRG ADV)/Personal_Estefania/Network/standard/toASint/")
class<-read.delim("~/Dropbox (CRG ADV)/Personal_Estefania/Network/summaryLinks/edgelist/class_colors_2020.txt", header = T)
##################################################################
expected<-read.table(
"../expected_alphSort.tab", sep="\t", header = T, row.names = 1)
checks<-which(is.na(match(expected$source, class$gene.name.VT))); length(checks)
checkt<-which(is.na(match(expected$target, class$gene.name.VT)));length(checkt)
remove<-union(checks, checkt)
length(remove)#786
#nodos que no estan en las redes pero estan en las DBs pues son interesantes
expected<-expected[-remove,]
dim(expected)
length(unique(expected$l11))
class(expected)
#load libraries
##############################
sourcedf<-as.data.frame(table(expected$short))
colnames(sourcedf)<-c("source","times")
head(sourcedf)
#figure 1
##################################################################
colnames(expected)
listSources<-list()
nodesSourceList<-list()
expectedTop<-expected
sources<-names(table(expectedTop$short))
length(unique(expectedTop$l11))
length(unique(expected$l11))
sources

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
nodesSourceList
#FIG 2: check Nodes in all the networks
#FIG 2b: check Nodes in all the networks
##################################################################
#evidence non from RNA-Seq
names(listSources)
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
setwd("../diffRho/ES_all/")
##################################################
edgelists<-read.table("all.txt", header = F, stringsAsFactors = F ,sep="\t")
#agrego las top 500:
finalDF<-data.frame()
linksList<-list()
length(linksList)
nodesList<-list()
dim(edgelists)#6

#All de caada network:
for (i in 1:nrow(edgelists)){
  file <- edgelists$V1[i]
  print(file)
  name <- gsub(".tab", "", edgelists$V1[i])
  print(name)
  edgeListLC <-  read.table(edgelists$V1[i], header = T, row.names = 1)
# edgeListLC <-  read.table(edgelists$V1[i], header = T, row.names = 1, nrows = 500)
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
###############################################
names(RNASeq[4])
linksListAll<-c(RNASeq[4], linksList)
names(linksListAll)
lapply(linksListAll, length)
getwd()
png(file="checkLinks_ESAll_vsLABCHIP.png",  width = 12, height = 8, units = 'in', res = 300)

upset(fromList(linksListAll),
      sets = c("Labchip", 
               "ES_all_dscaled_015",
                "ES_all_sscaled_02"), #500
      keep.order = TRUE,
      nsets = 3, nintersects = 10,
      order.by = "freq",
      text.scale = c(1,#eje Y
                     2,#eje y labels
                     2,#Set size
                     2,#set size ñabels
                     2,#dsrt lables
                     4#barplot labels
      ))
dev.off()
lapply(linksListAll, length)
names(linksListAll)
upset(fromList(linksListAll),
      sets = c("Labchip",  "ES_Q3_2_DS_st"))

upset(fromList(linksListAll),
      sets = c("Labchip",  "ES_Q3_3_DS_st"))

upset(fromList(linksListAll),
      sets = c("Labchip",  "ES_Q3_4_DS_st"))

upset(fromList(linksListAll),
      sets = c("Labchip",  "ES_Q3_5_DS_st"))

146/500*100
#################################################################

#############################################################
#cual de las networks recupera mejor los datos experimentales?
names(NoRNASeq)
experimental<-NoRNASeq
experimental[[4]]<-NULL
allExp<-list(unique(unlist(experimental)))
length(allExp[[1]])

rnaseq_exp<-c(linksListAll, allExp)
length(rnaseq_exp)#10
names(rnaseq_exp)[10]<-"allExp"
names(rnaseq_exp)
lapply(rnaseq_exp, length)

png(file="rnaseq_vs_exp_top500.png",  width = 8, height =4, units = 'in', res = 300)
  upset(fromList(rnaseq_exp),
        sets = c("Labchip", 
                 "ES_Q3_3_DS_st", #500
                 "ES_Q3_4_DS_st",#500
                 "ES_Q3_2_DS_st",#500
                 "ES_Q3_5_DS_st",
                 "allExp"), #500
       keep.order = TRUE,
        nsets = 6, nintersects = 10,
        order.by = "freq",
        text.scale = c(1,#eje Y
                       1,#eje y labels
                       1,#Set size
                       1,#set size ñabels
                       1,#dsrt lables
                       1#barplot labels
        ))
  dev.off()
#############################################################
#de acuerdo a esto, numeros:
lapply(rnaseq_exp, length)
upset(fromList(rnaseq_exp), sets = c("Labchip",  "allExp"))#237/541*100
upset(fromList(rnaseq_exp), sets = c( "ES_Q3_2_DS_st" , "allExp"))#210
upset(fromList(rnaseq_exp), sets = c( "ES_Q3_3_DS_st" , "allExp"))#197
upset(fromList(rnaseq_exp), sets = c( "ES_Q3_4_DS_st" , "allExp"))#215
upset(fromList(rnaseq_exp), sets = c( "ES_Q3_5_DS_st" , "allExp"))#214
#########################################
