library(ggplot2)
library(ggpubr)
#################################################################
setwd("SpliceNet/")
source("subsetFromGCDistribution.R")
workingDir<-"../SpliceNetData/data/"
outputPath<-"../SpliceNetRes/GC/"
## all posibble exons, changing: 16722
#this table was already prepared
inputFile<-"dPSI_EXONS.tab"
name="EXONS"
setwd(workingDir)
#############
table <-
  read.table(
    inputFile,
    sep = "\t",
    dec = "," ,
    header = T,
    row.names = 1
  )
print(paste("Total input", nrow(table), sep = ": "))
#########################################
features <- read.table("Exons_gc.tab", sep = "\t", header = T)
rownames(features) <- features$EVENT
#########################################
GCi <- match(rownames(table), rownames(features))
GC_Df <- data.frame(EVENT = row.names(features)[GCi],
                    GC = features$EXON_GCC[GCi])
table$GC <- features$EXON_GCC[GCi]
length(which(is.na(table$GC)))
table <- table[!is.na(table$GC), ]#53 NAs
###############################################################################
subsetFromGCDistribution(dPSIFile = table,
                         all = exons_gc,
                         name, outputPath)
###############################################################################
#forPlot:
exons_gc$changing <- rep("NC", nrow(exons_gc))
ee <- match(rownames(table), rownames(exons_gc))
changing <- exons_gc[ee, ]
dim(changing)#16668
dim(exons_gc)# 223487
changing$changing <- rep("CH", nrow(changing))

forPlot <- rbind(exons_gc, changing)

fileBP <- paste(outputPath,
                paste(name, "paired_boxplot.png", sep = "."), sep = "/")
png(
  fileBP,
  width = 12,
  height = 8,
  units = 'in',
  res = 300
)

ggplot(forPlot,  aes(x = changing, y = EXON_GCC)) +
  geom_boxplot(width = 0.1) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1,
    color = "red"
  ) +
  stat_compare_means(label.y = 1, label.x = 2) +
  theme_minimal()
dev.off()
##################################
#Distributions

fileDens <- paste(outputPath,
                  paste(name, "gc_density.png", sep = "."), sep = "/")
png(
  fileDens,
  width = 12,
  height = 8,
  units = 'in',
  res = 300
)

ggplot(forPlot,  aes(x = EXON_GCC)) +
  geom_density(aes(color = changing)) +
  theme_minimal()
dev.off()
