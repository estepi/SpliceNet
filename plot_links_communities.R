###############################################
library(igraph)
###############################################
folder<-"~/Dropbox (CRG ADV)//Personal_Estefania/SpliceNetRes/cor03/"
setwd(folder)

class<-read.delim("~/Dropbox (CRG ADV)//Personal_Estefania/Network/standard/diffRho/plot/class_colors_2020.txt", header = T)
head(class)

#############################################
ES<-read.table("ES_edgelist.tab", header = T); head(ES)
IR<-read.table("IR_edgelist.tab", header = T); head(IR)
A5<-read.table("A5_edgelist.tab", header = T); head(A5)
A3<-read.table("A3_edgelist.tab", header = T); head(A3)
###########################
g1 <- graph_from_data_frame(A3[,7:8], directed = FALSE)
#g1
length(E(g1))#
length(V(g1))#

layg1<-layout.fruchterman.reingold(g1)
rownames(layg1)<-V(g1)$name

cfg <- cluster_fast_greedy(g1,
                           modularity = TRUE,
                           weights = E(g1)$weight)

clp <- cluster_label_prop(g1, weights = E(g1)$weight)

colorColsLinks<-rep("lightsalmon", length(table(cfg$membership)))
colorColsLinks<-rep(NULL, length(table(cfg$membership)))
###classify links by communities
###############################################################
ii <- match(names(V(g1)), class$gene.name.VT)
V(g1)$label_color = class$color[ii]
V(g1)$label_color[is.na(V(g1)$label_color)] <- "grey"
###############################################################
dset = "A3_all"
pdfile <- paste(dset, "network.pdf", sep = "_")

pdf(pdfile, width = 12,  height = 12)
plot(
    cfg,
    g1,
    layout=layg1, 
    vertex.size = 4,
    edge.width = 1,
    edge.color = "grey",
    vertex.label.cex = 1,
    vertex.label.color = "black",
    vertex.label.dist = 0.8,
    main = dset,
    mark.col = colorColsLinks[-1],
    mark.border = colorColsLinks[-1],
    col = V(g1)$label_color)
dev.off()

####################################################################
