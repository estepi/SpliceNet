library(optparse)
library(utils)
library(igraph)
library(scales)
library(Hmisc)
library(tibble)
library(tidyr)
library(utils)
library(dplyr)
#####################################################################
option_list = list(
  make_option(
    c("-f", "--file"),
    type = "character",
    default = NULL,
    help = "eventscaled file name",
    metavar = "character"
  ),
  
  make_option(
    c("-m", "--min"),
    type = "numeric",
    default = "0.1",
    help = "min cor value"
  ),
  
  make_option(
    c("-p", "--fdr"),
    type = "numeric",
    default = "1",
    help = "max fdr value"
  ),
  
  make_option(
    c("-n", "--name"), 
    type = "character",
    default = "test",
    help = "file name"
    ),
  make_option(
    c("-b", "--bin"),
    type = "character",
    default = getwd(),
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

#####################################################################
print(opt$file)
name<-opt$name
print(paste("name:",name))
minCor<-opt$min
print(paste("Min cor:",minCor))
fdr<-opt$fdr
scripts<-opt$bin
sampleData<-read.table(opt$file, sep="\t", header=T) 
#####################################################################
# to run interactively
#t_Zvals<-read.table("")
#sampleData <- as.matrix(t_Zvals)
#sampleClass<-read.csv("class_colors.tab", sep="\t", header = T)
#name<-"ES"
#minCor=0.3
#fdr=0.1
# scripts<-".../SpliceNet"
#####################################################################
sc3<-paste(scripts, "functions_fdr_MC_cor.R", sep="/")
source(sc3)
M <- as.matrix(sampleData)
corM <- getCorM(M, min=0.4, fdr=1)
head(corM)
write.table(
  corM,
  file = paste(name, "edgelist.tab" , sep = "_"),
  sep = "\t",
  quote = F
)
#####################################################################
