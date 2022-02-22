#first generate data: TRUE and RANDOM DATA
library(optparse)# yes
library(parallel)
library(Hmisc)
library(tibble)
library(tidyr)
library(utils)
library(dplyr)
##############################################################
library(igraph)
library(scales)
library(reshape) 
library(ggplot2)
##############################################################
option_list = list(
  make_option(
    c("-f", "--file"),
    type = "character",
    default = NULL,
    help = "matrix file name",
    metavar = "character"
  ),
  make_option(
    c("-r", "--random"),
    type = "numeric",
    default = "10",
    help = "number of random matrix to test [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-s", "--start"),
    type = "numeric",
    default = "0.3",
    help = "value of rho to start scanning [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-e", "--end"),
    type = "numeric",
    default = "0.8",
    help = "value of rho to finish scanning [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-i", "--interval"),
    type = "numeric",
    default = "0.05",
    help = "intervale between each rho computation [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cores"),
    type = "numeric",
    default = "1",
    help = "number of cores to use[default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-b", "--bin"),
    type = "character",
    default = getwd(),
    help = "abs path folder of scripts [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-n", "--name"),
    type = "character",
    default = "test",
    help = "abs path folder of scripts [default= %default]",
    metavar = "character"
  )
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
##############################################
t_Zvals<-read.table(opt$file)
sampleData <- as.matrix(t_Zvals)
name<-opt$name
start<-opt$start
end<-opt$end
interval<-opt$interval
ncores <-opt$cores
fdr<-0.01
scripts<-opt$bin
##############################################
t_Zvals<-read.table("/home/estepi/Documents/SpliceNetData/input/ES_subset3500_sscaled.tab")
sampleData <- as.matrix(t_Zvals)
name<-"ES"
start<- 0.1
end<- 0.4
interval<-0.02
ncores <- 1
fdr<-0.01
scripts<-"/home/estepi/Documents/SpliceNet"

sc3<-paste(scripts, "functions_centrality_MC_cor.R", sep="/")
#esta funcion se carga de sc3
source(sc3)
start
end
DGList<-getCentralityByCOR(inputM, start, end, interval, fdr)
length(DGList)
#order according sample calssification

x[is.na(x)] <- 0
#hacer el heatmap
DGdf<-DGList[[1]]
DGdf[is.na[DGdf]]<-0

write.table(DGdf, "degree.tab", sep="\t", col.names = NA)

NDGdf<-DGList[[2]]
NDGdf[is.na[NDGdf]]<-0

write.table(NDGdf, "norm_degree.tab", sep="\t", col.names = NA)

sampleClass<-read.csv("~/Documents/summaryLinks/class_colors.tab", sep="\t", header = T)
head(sampleClass)

forPlot<-NDGdf
head(NDGdf)
forPlot<-DGdf

forPlot$gene<-rownames(forPlot)

forPlotMelt <- melt(forPlot, variable_name = "gene")
colnames(forPlotMelt)[2]<-"threshold"

ii<-match(forPlot$gene, sampleClass$gene.name.VT)
forPlotMelt$order<-sampleClass$ORDER[ii]
head(forPlotMelt)
forPlotMelt
forPlotMelt$rescale<-rescale(forPlotMelt$value, to=c(0.1,1))
ii<-match(forPlotMelt$gene, sampleClass$gene.name.VT)
forPlotMelt$class<-sampleClass$Class...family.small[ii]
forPlotMelt$color<-sampleClass$color[ii]
head(forPlotMelt)
colnames(forPlotMelt)


complex<-unique(forPlotMelt$class)
head(complex)
cc<-match(complex, sampleClass$Class...family.small)
color<-sampleClass$color[cc]
names(color)<-complex
color
head(forPlotMelt)

#################################################
pdf("Degree.pdf", width=6, height=12)
ggplot(data=forPlotMelt) +
  geom_tile(aes(x=threshold, 
                y=reorder(gene,order),
                fill=class, 
                alpha=rescale)) + 
  theme_classic()+
  scale_fill_manual(values=color) + #assign tissues colors
  scale_alpha_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))+
  theme(axis.text.y = element_text( size=3))+
  theme(legend.position = "none") 
dev.off()
#################################################






