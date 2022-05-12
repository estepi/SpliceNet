#first generate data: TRUE and RANDOM DATA
library(igraph)
library(scales)
###############################################################################
NetworkDesc <- function(g, name)
{
  ###################################################
  DG <- degree(g)
  NDG <- DG / length(E(g))#
  PR <- page.rank(g)$vector  		#Pagerank Score
  OPR <- order(PR, decreasing = TRUE)
  BC <- betweenness(g)		#Betweeness Centrality
  OBC <- order(BC, decreasing = TRUE)
  CC <- closeness(g)		#Closeness Centrality
  OCC <- order(CC, decreasing = TRUE)
  RANK <- rank(PR) + rank(BC) + rank(CC)
  ORANK <-
    order(RANK, decreasing = TRUE)	#Ordering of Nodes based on aggregate score.
  file6 <- paste(name, "CentralityRanking.tab", sep = "_")
  df3 <- data.frame(DG, NDG, PR, OPR, BC, OBC, CC, OCC, RANK, ORANK)
  head(df3)
  write.table(df3,
              file6,
              col.names = NA,
              quote = F,
              sep = "\t")
}
#####################################################################
setwd("SpliceNetRes/GC/networks/cor03/")

file.names <- dir(path = ".", pattern = "\\_edgelist.tab")
file.names

for (i in 1:length(file.names))
{
  name <- gsub("_edgelist.tab", "", file.names[i])
  print(name)
  head(read.table(file.names[i]))
  g <- graph_from_data_frame(read.table(file.names[i])[, 7:8],
                             directed = FALSE)
  
  NetworkDesc(g, name)
}



      


