###############################################
library(igraph)
###############################################
folder <- "SpliceNetRes"
setwd(folder)
class <-
    read.delim("../SpliceNetData/class_colors_2020.txt", header = T)
head(class)
#############################################
Q1 <- read.table("Q1_edgelist.tab", header = T)
head(Q1)
dim(Q1)#453
Q1top <- Q1[Q1$absCor > 0.4, ]

g1 <- graph_from_data_frame(Q1top[, 7:8], directed = FALSE)
onlyG1 <- difference(g1, g3)
#g1

layg1 <- layout.fruchterman.reingold(onlyG1)
rownames(layg1) <- V(onlyG1)$name

cfg <- cluster_fast_greedy(onlyG1,
                           modularity = TRUE,
                           weights = E(onlyG1)$weight)

clp <- cluster_label_prop(onlyG1, weights = E(onlyG1)$weight)

colorColsLinks <- rep("lightsalmon", length(table(cfg$membership)))
colorColsLinks <- rep(NULL, length(table(cfg$membership)))
###classify links by communities
###############################################################
ii <- match(names(V(onlyG1)), class$gene.name.VT)
V(onlyG1)$label_color = class$color[ii]
V(onlyG1)$label_color[is.na(V(onlyG1)$label_color)] <- "grey"
###############################################################
dset = "onlyG1"
pdfile <- paste(dset, "network.pdf", sep = "_")

pdf(pdfile, width = 12,  height = 12)
plot(
    cfg,
    onlyG1,
    layout = layg1,
    vertex.size = 4,
    edge.width = 1,
    edge.color = "grey",
    vertex.label.cex = 1,
    vertex.label.color = "black",
    vertex.label.dist = 0.8,
    main = dset,
    mark.col = colorColsLinks[-1],
    mark.border = colorColsLinks[-1],
    col = V(onlyG1)$label_color
)
dev.off()
getwd()
##########################################33
Q3 <- read.table("Q3_edgelist.tab", header = T)
head(Q3)
dim(Q3)#453
Q3top <- Q3[Q3$absCor > 0.4, ]
dim(Q3top)#481
###########################
g3 <- graph_from_data_frame(Q3top[, 7:8], directed = FALSE)
onlyG3 <- difference(g3, g1)
length(E(onlyG3))#
length(V(onlyG3))#

layg3 <- layout.fruchterman.reingold(onlyG3)
rownames(layg3) <- V(onlyG3)$name

cfg <- cluster_fast_greedy(onlyG3,
                           modularity = TRUE,
                           weights = E(onlyG3)$weight)
clp <- cluster_label_prop(g3, weights = E(onlyG3)$weight)

colorColsLinks <- rep("lightsalmon", length(table(cfg$membership)))
colorColsLinks <- rep(NULL, length(table(cfg$membership)))
###############################################################
ii <- match(names(V(onlyG3)), class$gene.name.VT)
V(onlyG3)$label_color = class$color[ii]
V(onlyG3)$label_color[is.na(V(onlyG3)$label_color)] <- "grey"
###############################################################
dset = "onlyG3"
pdfile <- paste(dset, "network.pdf", sep = "_")

pdf(pdfile, width = 12,  height = 12)
plot(
    cfg,
    onlyG3,
    layout = layg3,
    vertex.size = 4,
    edge.width = 1,
    edge.color = "grey",
    vertex.label.cex = 1,
    vertex.label.color = "black",
    vertex.label.dist = 0.8,
    main = dset,
    mark.col = colorColsLinks[-1],
    mark.border = colorColsLinks[-1],
    col = V(onlyG3)$label_color
)
dev.off()

