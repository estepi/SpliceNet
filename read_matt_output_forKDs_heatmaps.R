library(corrplot)
library(pheatmap)
library(igraph)
library(RColorBrewer)
###############################################
output <- "SpliceNetRes/"
#define where you want to put output
dir(output)
folder <- "SpliceNetRes"
#in folder there are a folder for each condition and inside should be the file:
#overview_of_features_and_comparisons.tab
pname <- "Exons_10_bgC" #is the name of the project
pattern <- "exons_do_c" #folders to be analyzed
outfile <- paste(output, name, sep = "/")
setwd(folder)
kds <-
  data.frame(dir(folder, pattern = pattern), stringsAsFactors = F)
colnames(kds) <- "kd"
featuresTable <- data.frame()

for (i in 1:nrow(kds))
{
  dset <- kds$kd[i]
  print(dset)
  name <- gsub("exons_do_c_", "", kds$kd[i])#change for project name
  print(name)
  featuresTable <- rbind(featuresTable,
                         data.frame(t(
                           read.table(
                             file = paste(dset, "overview_of_features_and_comparisons.tab", sep = "/"),
                             header = T,
                             row.names = 1,
                             sep = "\t"
                           )
                         ))
                         ["PVAL_MANN_WHITNEY_U_TEST_bgConst_VS_do", ])
  head(featuresTable)
  rownames(featuresTable)[i] <- name
}
###########################################################
head(featuresTable)  
colnames(featuresTable)
rownames(featuresTable)
###########################################################
Numeric <-
  data.frame(matrix(
    unlist(apply(featuresTable, 1, as.numeric)),
    ncol = ncol(featuresTable),
    nrow = nrow(featuresTable),
    byrow = T
  ))
colnames(Numeric) <- colnames(featuresTable)
rownames(Numeric) <- rownames(featuresTable)
corKDs <- cor(t(Numeric))
corFeatures <- cor(Numeric)
head(Numeric)
###############################################
#agregamos color por communities:
links <- read.table(
  "edgelist.tab",
  stringsAsFactors = F,
  row.names = 1,
  header = T
)
rownames(links) <- paste(links$source, links$target, sep = "-")
head(links)
g <- graph_from_data_frame(links[, 7:8], directed = FALSE)
net <- g
###############################################
sampleClass <-
  read.delim("SpliceNetData/class_colors_2020.txt",
             header = T)
ii <- match(V(net)$name, sampleClass$gene.name.VT)
V(net)$color <- sampleClass$color[ii]
##################################################################
layout_with_dh <- layout.fruchterman.reingold(net)
cfg <-
  cluster_fast_greedy(net, modularity = TRUE, weights = E(net)$weight)
clp <- cluster_label_prop(net, weights = E(net)$weight)
###classify links by communities
nodes <- data.frame(node = V(net)$name)
rownames(nodes) <- nodes$node
iis <- match(nodes$node, cfg$names)
nodes$com <- cfg$membership[iis]
table(cfg$membership)#10 commu

set.seed(1234)
dset="ES"
######################################################
#subset KDs que tienen pval:
colsii <- match(rownames(nodes), colnames(corKDs))
rownames(nodes)[is.na(colsii)]

rowsii <- match(rownames(nodes), rownames(corKDs))
rownames(nodes)[is.na(rowsii)]

colsii <- colsii[!is.na(colsii)]
rowsii <- rowsii[!is.na(rowsii)]

corKDs_Com <- corKDs[rowsii, colsii]
####################################
paletteLength = 10
myColor <-
  colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(
  seq(min(corKDs_Com), 0,
      length.out = ceiling(paletteLength / 2) + 1),
  seq(
    max(corKDs_Com) / paletteLength,
    max(corKDs_Com),
    length.out = floor(paletteLength / 2)
  )
)
#############################################
newCols <-
  colorRampPalette(grDevices::rainbow(length(unique(nodes$com))))
mycolors <- newCols(length(unique(nodes$com)))
names(mycolors) <- unique(nodes$com)
mycolorsL <- list(com = mycolors)
mycolors
#############################################
pdfile <- paste(dset, "network.pdf", sep = "_")
pdf(pdfile, width = 10,  height = 10)
plot(
  cfg,
  net,
  vertex.size = 10,
  edge.width = 5,
  edge.color = "blue",
  vertex.label.cex = 0,
  vertex.label.color = "black",
  layout = layout_with_dh * 0.5,
  main = dset,
  mark.col = mycolors,
  mark.border = mycolors,
  col = V(net)$color
)
dev.off()
############################################
nodes$node <- NULL
head(corKDs_Com)
dim(corKDs_Com)

pdfile <- paste(dset, "heatmap.pdf", sep = "_")
pdf(pdfile, width = 10,  height = 10)
pheatmap(
  corKDs_Com,
  fontsize_col = 10,
  fontsize_row = 10,
  color = myColor,
  breaks = myBreaks,
  annotation_row = nodes,
  annotation_col = nodes,
  annotation_colors = mycolorsL
)
dev.off()

