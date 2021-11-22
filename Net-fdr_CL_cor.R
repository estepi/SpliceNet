#first generate data: TRUE and RANDOM DATA
library(optparse)# yes
library(parallel)
library(Hmisc)
library(tibble)
library(tidyr)
library(utils)
library(dplyr)
##############################################################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="matrix file name", metavar="character"),
  make_option(c("-r", "--random"), type="numeric", default="10", 
              help="number of random matrix to test [default= %default]", metavar="character"),
  make_option(c("-s", "--start"), type="numeric", default="0.3", 
              help="value of rho to start scanning [default= %default]", metavar="character"),
  make_option(c("-e", "--end"), type="numeric", default="0.8", 
              help="value of rho to finish scanning [default= %default]", metavar="character"),
  make_option(c("-i", "--interval"), type="numeric", default="0.05", 
              help="intervale between each rho computation [default= %default]", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default="1", 
              help="number of cores to use[default= %default]", metavar="character"),
  make_option(c("-b", "--bin"), type="character", default=getwd(), 
              help="abs path folder of scripts [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="test",
             help="abs path folder of scripts [default= %default]", metavar="character")
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
NumRandomM<-opt$random
start<-opt$start
end<-opt$end
interval<-opt$interval
ncores <-opt$cores
fdr<-0.01

print(c(start,end, interval, NumRandomM))
scripts<-opt$bin

sc3<-paste(scripts, "functions_fdr_MC_cor.R", sep="/")
#esta funcion se carga de sc3
source(sc3)

estimateFDR(  inputM=sampleData,  
              NumRandomM,  start,  end,  interval, ncores, name, fdr)


