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
sc3<-paste(scripts, "functions_fdr_MC_cor.R", sep="/")
#esta funcion se carga de sc3
source(sc3)


sampleData<-read.table(opt$file, sep="\t", header=T) 
M <- as.matrix(sampleData)
corM <- getCorM(M, minCor, fdr)
write.table(corM, 
	   file=paste(name, "edgelist.tab" , sep="_"), 
	   sep="\t",quote=F)

