##############################################################
library(parallel)##multicores
library(Hmisc)# correlation matrix
library(optparse)#parse input 
library(tibble)#rownmaes_to_columns
library(tidyr)#gather
library(dplyr)#left join
library(igraph)# networks
library(scales)#scaling
library(reshape) #for plot
library(ggplot2)#for plot
#


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
    c("-o", "--order"),
    type = "character",
    default = NULL,
    help = "classification file name",
    metavar = "character"
  ),
  make_option(
    c("-s", "--start"),
    type = "numeric",
    default = "0.1",
    help = "value of rho to start scanning [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-e", "--end"),
    type = "numeric",
    default = "0.4",
    help = "value of rho to finish scanning [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-i", "--interval"),
    type = "numeric",
    default = "0.02",
    help = "interval between each cor value [default= %default]",
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
    help = "sufix for output [default= %test]",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
##############################################
print(opt$file)
print(opt$order)
print(opt$name)
print(opt$start)
print(opt$end)
print(opt$interval)
print(opt$cores)
print(opt$bin)
##############################################
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
##############################################
t_Zvals<-read.table(opt$file)
sampleData <- as.matrix(t_Zvals)
sampleClass<-read.delim(opt$order, sep="\t", header = T)
print(head(sampleClass))
print(colnames(sampleClass))
print(head(sampleClass$color))
name<-opt$name
start<-opt$start
end<-opt$end
interval<-opt$interval
ncores <-opt$cores
fdr<-0.01
scripts<-opt$bin
sc3 <- paste(scripts, "functions_centrality_MC_cor.R", sep = "/")
#esta funcion se carga de sc3
source(sc3)
DGList <- getCentralityByCOR(sampleData, start, end, interval, fdr)
#########################################################
print("Finish computation. Lets plot")
FileName<-paste(name, start, end, sep="-")
DGdf <- DGList[[1]]
DGdf[is.na(DGdf)] <- 0
degreeFile <- paste(FileName, "degree.tab", sep = "_")
write.table(DGdf, degreeFile, sep = "\t", col.names = NA)
#########################################################
NormDegreeFile <- paste(FileName, "norm_degree.tab", sep = "_")
NDGdf <- DGList[[2]]
NDGdf[is.na(NDGdf)] <- 0
write.table(NDGdf, NormDegreeFile, sep = "\t", col.names = NA)
#########################################################
