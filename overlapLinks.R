library(ggplot2)
library(UpSetR)
library(igraph)
library(ggpubr)
##################################################################
setwd("~/Documents/SpliceNetRes/cor03/")
class<-read.delim("~/Documents/summaryLinks/Q3_R_SS_Networks/class_colors.txt", header = T)
##################################################################
expected<-read.table(
"/home/estepi/Documents/summaryLinks/Q3_R_SS_Networks/overlap_vs_expected/expected_alphSort.tab", sep="\t", header = T, row.names = 1)
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

for ( i in 1:length(sources))
{
  print(sources)[i]
  listSources[[i]]  <-expected$l11[expected$short==sources[i]]
  names(listSources)[i]<-sources[i]
  #######################################
  n1<-unique(c(as.character(expected$source[expected$short==sources[i]]),
               as.character(expected$target[expected$short==sources[i]])))
  nodesSourceList[[i]]<-n1
  names(nodesSourceList)[[i]]<-sources[i]
  
}

nodesSourceList
#FIG 2: check Nodes in all the networks
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
labchiplinks<-unlist(RNASeq[4])
length(labchiplinks)#541
length(unique(onlyRNASeqLinks))#10291
length(which(labchiplinks%in%onlyRNASeqLinks))
##################################################
#convierto en la lista de redes:
edgelists<-read.table("edgelist.tab", header = F, stringsAsFactors = F ,sep="\t")
edgelists
edgelists[1:3,]; 
dim(edgelists)#6
finalDF<-data.frame()
linksList<-list()
length(linksList)
nodesList<-list()
edgelists#5

for (i in 1:nrow(edgelists))
{ 

    file<-edgelists$V1[i]
    print(file)
    name<-gsub("_edgelist.tab", "", edgelists$V1[i])
    print(name)
    edgeListLC<-read.table(edgelists$V1[i], header = T, row.names = 1);
    print(head(edgeListLC));
    print(dim(edgeListLC))
    ###########################################################
    #Links hay que hacer el sort
    edgeListLCSort<-data.frame(edgeListLC$source, edgeListLC$target)
    colnames(edgeListLCSort)<-c("source", "target")
    l11<-paste(edgeListLCSort$source, edgeListLCSort$target, sep="-")
    linksList[[i]]<-l11
    names(linksList)[[i]]<-name
    ###########################################################
    #nodes
    n1<-unique(c(as.character(edgeListLC$source),
               as.character(edgeListLC$target)))
    nodesList[[i]]<-n1
    names(nodesList)[[i]]<-name
}
##################################################################
lapply(linksList, head)
allGenes<-unlist
names(nodesList)
df<-data.frame(names=unlist(lapply(nodesList, length)))
head(df)
names(nodesList)
length(nodesList)#6, saco el ALL y agrego total /singel

nodesListSinCore<-nodesList
nodesListSinCore[1]<-NULL
length(nodesListSinCore)#5
#FIG5 Check NODES
###################################
#chequear que KDs son nuevos y cuáles son labchip
head(class)
dim(class)#305
class$Gene.ID
class$Gene.Symbol
###################################
names(linksList)
ESNet<-linksList
ESNet[1]<-NULL
names(ESNet)
names(NoRNASeq)
NoRNASeq
checkLinks<-c(linksList[6],NoRNASeq[4])
lapply(checkLinks, length)
#fig6
png(file="checkKDsNames_5E.png",  width = 12, height = 8, units = 'in', res = 300)
upset(fromList(checkLinks), 
      order.by = "freq",
      text.scale = c(1,#eje Y
                     2,#eje y labels
                     2,#Set size
                     2,#set size ñabels
                     2,#dsrt lables
                     5#barplot labels
                     ))
dev.off()
###################################
names(linksList)
linksList[1]<-NULL
lapply(linksList, length)
#FIG7
png(file="checkLinks.png",  width = 12, height = 8, units = 'in', res = 300)
upset(fromList(linksList),
      sets = c("labC_ST", 
               "ES_DS_total",
               "ES_SS_total",
               "ES_DS_single",
               "ES_SS_single"),
      keep.order = TRUE,
      nsets = 5, nintersects = 15,
      order.by = "freq",
      text.scale = c(1,#eje Y
                     2,#eje y labels
                     2,#Set size
                     2,#set size ñabels
                     2,#dsrt lables
                     4#barplot labels
      ))
dev.off()
226/541
names(ESNet)
ESNet[1]
allLinksNets<-unique(unlist(ESNet))
neverRecovered<-labchiplinks[which(!labchiplinks%in%allLinksNets)]
length(neverRecovered)#226
#FIG7b
write(neverRecovered, "neverRecovered.tab")
SFs<-data.frame(table(unlist(strsplit(neverRecovered, "-"))))
head(SFs)

SFsTOP<-SFs[order(SFs$Freq,decreasing = T)[1:20],]
head(SFsTOP)
#FIG7c
png(width = 800, height = 400, file="neverRecovered_Vs_labchip.png")
ggplot(data=SFsTOP, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust=1))
dev.off()
##################################################################
neverRecoveredList<-list()
neverRecoveredList[[1]]<-neverRecovered
names(neverRecoveredList)<-"neverRecovered"
length(neverRecoveredList)
expected<-NoRNASeq[-4]
names(expected)
links<-c(expected, neverRecoveredList)
names(links)
length(links)#8
#FIG7d
png(file="NeverRecovered_upsetR.png",  width = 12, height =8, units = 'in', res = 300)
upset(fromList(links), order.by = "freq", 
      nintersects=20, nsets = 8,
      text.scale = c(1,#eje Y
                     2,#eje y labels
                     2,#Set size
                     2,#set size ñabels
                     2,#dsrt lables
                     2))
dev.off()
#################################################################
#all networks RNA-Seq vs experimental
#################################################################
all<-unique(unlist(links))
names(links)
auxdf<-data.frame(
  Experiental=all%in%links[[1]],
  GeneFusion=all%in%links[[2]],
  Homology=all%in%links[[3]],
  Locus=all%in%links[[4]],
  Phylo=all%in%links[[5]],
  Structure=all%in%links[[6]],
  TextMininh=all%in%links[[7]],
  NotRecovered=all%in%links[[8]])
rownames(auxdf)<-all
dim(auxdf)
#################################################################
#links experimental not recovered #41
colnames(auxdf)
expNotRecovered<-rownames(auxdf)[auxdf$Experiental==TRUE &
                       auxdf$NotRecovered ==TRUE]
length(expNotRecovered)#57 41+13
write(expNotRecovered, "expNotRecovered.tab")
#################################################################
evidence<-neverRecovered[neverRecovered%in%unlist(NoRNASeq)]
length(evidence)#225
head(evidence)
length(neverRecovered)#226
191+226
write(evidence, "neverRecoveredEvidence.tab")
#####################################################
#fig7e
SFsEvidence<-data.frame(table(unlist(strsplit(evidence, "-"))))
head(SFsEvidence)
SFsEvidenceTop<-SFsEvidence[order(SFsEvidence$Freq, decreasing = T)[1:20],]
head(SFsEvidenceTop)
#####################################################
#fig7e
png(width = 600, height = 400, file="neverRecoveredEvidenceTop.png")
ggplot(data=SFsEvidenceTop, aes(x=reorder(Var1, -Freq), y=Freq)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust=1))
dev.off()
#####################################################
####################################################
#extraer los non recoverd por labchip que estan en xperimentl
#NEW HC
#FIG8
all<-unique(unlist(linksList))
names(linksList)
length(all)
auxdf<-data.frame(ES_DS_single=all%in%linksList[[1]],
                  ES_DS_total=all%in%linksList[[2]],
                  ES_SS_single=all%in%linksList[[3]],
                  ES_SS_total=all%in%linksList[[4]],
                  labC=all%in%linksList[[5]])

rownames(auxdf)<-all
length(which(auxdf$labC==TRUE))#541
head(auxdf)
new<-rownames(auxdf)[auxdf$ES_DS_single==TRUE &
                  auxdf$ES_DS_total==TRUE&
                  auxdf$ES_SS_single==TRUE &
                    auxdf$ES_SS_total==TRUE &
                  auxdf$labC==FALSE ]
length(new)#488

#compare with experimental
newList<-list()
newList[[1]]<-new
names(newList)<-"new"
length(newList)
newlinks<-c(expected, newList)
names(newlinks)
names(expected)
length(newlinks)#8
#FIG8
png(file="new_upsetR.png",  width = 12, height =8, units = 'in', res = 300)
upset(fromList(newlinks), 
      order.by = "freq",  
      nintersects=20, 
      nsets = 8,
      text.scale = c(1,#eje Y
                     2,#eje y labels
                     2,#Set size
                     2,#set size ñabels
                     2,#dsrt lables
                     2))#barplot labels

dev.off()

write(new, "new.tab")
#Sfs?
notCo<-notLC[notLC%in%labchipNames$gene.name.VT]
length(notCo)#67
newSamples<-notLC[!notLC%in%labchipNames$gene.name.VT]
length(newSamples)

SFsNew<-data.frame(table(unlist(strsplit(new, "-"))))
table(SFsNew$class)
dim(SFsNew)#170

SFsNew$class<-rep("old")
SFsNew$class[match(SFsNew$Var1, notCo)]<-"LabchipNotCon"
SFsNew$class[match(SFsNew$Var1, newSamples)]<-"newSamples"

table(SFsNew$class)

SFsNewTop<-SFsNew[order(SFsNew$Freq, decreasing = T)[1:20],]
head(SFsNewTop)
unconnected
new


#####################################################
#FIG8b
png(width = 600, height = 400, file="neverRecoveredEvidenceTop.png")
ggplot(data=SFsNewTop, aes(x=reorder(Var1, -Freq), y=Freq, fill=class)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust=1))
dev.off()
#####################################################
#####################################################
#compare with RNA-Seq
names(RNASeq)

CrossRNASeq<-RNASeq[1:3]

newlinksCross<-c(CrossRNASeq, newList)
names(newlinksCross)
length(newlinks)#8
#FIG9
png(file="new_vs_cross.png",  width = 6, height =4, units = 'in', res = 300)
upset(fromList(newlinksCross), order.by = "freq",  nintersects=20, nsets = 4)
dev.off()
###################################################
all<-unique(unlist(newlinksCross))
length(all)
head(all)
names(newlinksCross)
auxdf<-data.frame(CROSS_AS=all%in%newlinksCross[[1]],
                  CROSS_GE_Down=all%in%newlinksCross[[2]],
                  CROSS_GE_Up=all%in%newlinksCross[[3]],
                  new=all%in%newlinksCross[[4]])
head(auxdf)
rownames(auxdf)<-all

head(auxdf)
newCrossAS<-rownames(auxdf)[auxdf$CROSS_AS==TRUE &
                       auxdf$CROSS_GE_Down==FALSE&
                       auxdf$CROSS_GE_Up==FALSE &
                       auxdf$new==TRUE]
length(newCrossAS)#87
write(newCrossAS, "newCrossAS.tab")

###############check in networks
#############################################################
#cual de las networks recupera mejor los datos experimentales?
names(linksList)
names(newlinks)
allExp<-list(unlist(newlinks))
head(allExp)
length(newlinks)
newlinks[[8]]<-NULL

rnaseq_exp<-c(linksList, allExp)
length(rnaseq_exp)
names(rnaseq_exp)[6]<-"allExp"
names(rnaseq_exp)

png(file="rnaseq_vs_exp.png",  width = 6, height =4, units = 'in', res = 300)
  upset(fromList(rnaseq_exp),
        sets = c("ES_DS_total",
                 "ES_SS_total",
                 "ES_DS_single",
                 "ES_SS_single",
                 "labC_ST",
                 "allExp"),
        keep.order = TRUE,
        nsets = 6, nintersects = 15,
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
# labC_ST [1] 541
upset(fromList(rnaseq_exp), sets = c("labC_ST",  "allExp"))
237/541*100
# ES_SS_single   [1] 8344
upset(fromList(rnaseq_exp), sets = c("ES_SS_single",  "allExp"))
459/8344*100
# ES_DS_single   [1] 2264
upset(fromList(rnaseq_exp), sets = c("ES_DS_single",  "allExp"))
303/2264 *100
# ES_SS_total   [1] 1036
upset(fromList(rnaseq_exp), sets = c("ES_SS_total",  "allExp"))
200/1036*100
# ES_DS_total   [1] 6021
upset(fromList(rnaseq_exp), sets = c("ES_DS_total",  "allExp"))
505/6021 * 100
# allExp  [1] 2596
upset(fromList(rnaseq_exp), sets = c("ES_DS_total",  "allExp"))
    
#domde estan los links en la frecuencia que tienen  validacion experimental?
#lo ultimo es ver, si modifico los thhreshold de las networks, para tener menos links y mas confiables
ES_DS_total<-read.table("../normalized/ES_DS_total.tab", row.names = 1)
head(ES_DS_total)
iid<-match(unlist(allExp), rownames(ES_DS_total))
length(which(is.na(iid)))#0
ES_DS_total$color<-rep("all", nrow(ES_DS_total))
ES_DS_total$color[iid]<-"exp"
colnames(ES_DS_total)
ES_DS_totalSort<-ES_DS_total[order(ES_DS_total$occur, decreasing = T),]
head(ES_DS_totalSort)
ES_DS_totalSort$Xvar<-1:nrow(ES_DS_totalSort)
head(ES_DS_totalSort)
table(ES_DS_totalSort$color)
###############################
png("ES_DS_total_occur_all.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_DS_totalSort, aes(Xvar, occur, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("total DS") 
dev.off()
###############################
png("ES_DS_total_occur_exp.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_DS_totalSort[ES_DS_totalSort$color=="exp",], aes(Xvar, occur, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("total DS exp") 
dev.off()
###############################
png("ES_DS_total_vs_exp_violin.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_DS_totalSort, aes(x=color, y=occur))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")    + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
    stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  ggtitle("total DS") 
dev.off()
###############################
ES_SS_t<-read.table("../normalized/ES_SS_total.tab", row.names = 1)
head(ES_SS_t)
iid<-match(unlist(allExp), rownames(ES_SS_t))
length(which(is.na(iid)))#0
ES_SS_t$color<-rep("all", nrow(ES_SS_t))
ES_SS_t$color[iid]<-"exp"
colnames(ES_SS_t)
ES_SS_tSort<-ES_SS_t[order(ES_SS_t$occur, decreasing = T),]
head(ES_SS_tSort)
ES_SS_tSort$Xvar<-1:nrow(ES_SS_tSort)
head(ES_SS_tSort)
table(ES_SS_tSort$color)
###############################
png("ES_SS_total_occur_all.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_tSort, aes(Xvar, occur, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("total SS") 
dev.off()
###############################
png("ES_SS_t_occur_exp.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_tSort[ES_SS_tSort$color=="exp",], aes(Xvar, occur, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("total SS exp") 
dev.off()
###############################
png("ES_SS_t_vs_exp_violin.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_tSort, aes(x=color, y=occur))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")    + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  ggtitle("total SS")
dev.off()
###############################
ES_SS_s<-read.table("../normalized/ES_SS_single.tab", row.names = 1)
head(ES_SS_s)
iid<-match(unlist(allExp), rownames(ES_SS_s))
length(which(is.na(iid)))#0
ES_SS_s$color<-rep("all", nrow(ES_SS_s))
ES_SS_s$color[iid]<-"exp"
colnames(ES_SS_s)
ES_SS_sSort<-ES_SS_s[order(ES_SS_s$abscor, decreasing = T),]
head(ES_SS_sSort)
ES_SS_sSort$Xvar<-1:nrow(ES_SS_sSort)
head(ES_SS_sSort)
table(ES_SS_sSort$color)
################################################
png("ES_SS_single_cor_all.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort, aes(Xvar, abscor, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("SS single") 
dev.off()
###############################
png("ES_SS_s_cor_exp.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort[ES_SS_sSort$color=="exp",], 
       aes(Xvar, abscor, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("SS single exp") 
dev.off()
###############################
png("ES_SS_s_vs_exp_violin.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort, aes(x=color, y=abscor))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")    + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  ggtitle("SS single")
dev.off()
###############################

###############################
ES_DS_s<-read.table("../normalized/ES_DS_single.tab")

head(ES_DS_s)
iid<-match(unlist(allExp), rownames(ES_DS_s))
length(which(is.na(iid)))#0
ES_DS_s$color<-rep("all", nrow(ES_DS_s))
ES_DS_s$color[iid]<-"exp"
colnames(ES_DS_s)
ES_DS_sSort<-ES_DS_s[order(ES_DS_s$abscor, decreasing = T),]
head(ES_DS_sSort)
ES_DS_sSort$Xvar<-1:nrow(ES_DS_sSort)
head(ES_DS_sSort)
table(ES_DS_sSort$color)
################################################
png("ES_SS_single_cor_all.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort, aes(Xvar, abscor, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("SS single") 
dev.off()
###############################
png("ES_SS_s_cor_exp.png",width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort[ES_SS_sSort$color=="exp",], 
       aes(Xvar, abscor, color=color))+
  geom_point()+
  theme_bw()+
  ggtitle("SS single exp") 
dev.off()
###############################
png("ES_SS_s_vs_exp_violin.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(ES_SS_sSort, aes(x=color, y=abscor))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")    + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  ggtitle("SS single")
dev.off()
###############################


###############################
#new + crossAS:
allExp
iid<-match(allExp, rownames(ES_DS_total))
length(which(is.na(iid)))#0

ES_DS_total$color<-rep("all", nrow(ES_DS_total))
ES_DS_total$color[iid]<-"exp"
colnames(ES_DS_total)
ES_DS_totalSort<-ES_DS_total[order(ES_DS_total$occur, decreasing = T),]
head(ES_DS_totalSort)
ES_DS_totalSort$Xvar<-1:nrow(ES_DS_totalSort)
head(ES_DS_totalSort)
table(ES_DS_totalSort$color)

ggplot(ES_DS_totalSort[ES_DS_totalSort$color=="newCrossAS",], aes(Xvar, occur, color=color))+
  geom_point()+
  theme_bw()
