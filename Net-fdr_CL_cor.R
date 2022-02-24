library(optparse)# yes
###############################################3
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
  stop("Example 
       Rscript -f ES_subset3500_sscaled.tab -r 2 -s 0.1 -e 0.2 -i 0.05 -c 1 -b ~/Documents/SpliceNet -n testCL 
       At least one argument must be supplied (input file).n", call.=FALSE)
}
##############################################
library(parallel)
library(Hmisc)
library(tibble)
library(tidyr)
library(utils)
library(dplyr)

t_Zvals<-read.table(opt$file)
sampleData <- as.matrix(t_Zvals)
name<-opt$name
NumRandomM<-opt$random
start<-opt$start
end<-opt$end
interval<-opt$interval
ncores <-opt$cores
fdr<-0.01
##############################################
#to run interactively, uncomment following lines:
#t_Zvals<-read.table("/home/estepi/Documents/SpliceNetData/input/ES_subset3500_sscaled.tab")
#sampleData <- as.matrix(t_Zvals)
#name<-"ES"
#start<- 0.1
#end<- 0.2
#interval<-0.02
#ncores <- 1
#fdr<-0.01
#NumRandomM<-5
#fdr<-0.01
##############################################
scripts<-"/home/estepi/Documents/SpliceNet"
sc3<-paste(scripts, "functions_fdr_MC_cor.R", sep="/")
source(sc3)

print(c(start,end, interval, NumRandomM))

estimateFDR(  inputM=sampleData,  
              NumRandomM,   start,  end,   interval,  ncores,  name, fdr)


