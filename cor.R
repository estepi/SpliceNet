library(optparse)
library(utils)
library(QUIC)
library(igraph)
library("RColorBrewer")
library(scales)
source("/no_backup/jvalcarcel/emancini/Network/scripts/CRobCor.R")
source("/no_backup/jvalcarcel/emancini/Network/scripts/functions_fdr_noMC.R")
source("/no_backup/jvalcarcel/emancini/Network/scripts/CreateGraph.R")
source("/no_backup/jvalcarcel/emancini/Network/scripts/CentralityRanking.R")
source("/no_backup/jvalcarcel/emancini/Network/scripts/Vscale.R")
##############################################################################
NetworkDesc<-function(g, name, Cdouble)
  {
df2<-data.frame(get.edgelist(g))
Cdouble<-as.data.frame(Cdouble)

print(paste("NumOfV:",length(V(g))), sep=":")
print(paste("NumOfE:",length(E(g))), sep=":")

edgeListLCSort<-data.frame(t(unlist(apply(df2,1, function(x){sort(x)}))))
head(edgeListLCSort)

final<-data.frame("source"=  edgeListLCSort$X1, 
                  "target" =edgeListLCSort$X2)
cor<-c()

for (i in 1:nrow(final)) {

ss<-as.character(final$source[i])
print(ss)
tt<-as.character(final$target[i])
print(tt)
cor[i]<-  Cdouble[ss,tt]
print(cor[i])
}

final$cor<-cor
final$abscor<-abs(cor)
file5<-paste(name, "edgelist.tab", sep="_")
write.table(final, file5, col.names = NA, quote = F, sep="\t")

#file6<-paste(name, "matrix.tab", sep="_")
#write.table(Cdouble, file6, col.names = NA, quote = F, sep="\t")

}

#####################################################################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="eventscaled file name", metavar="character"),
  make_option(c("-r", "--rho"), type="numeric", default="0.1", 
              help="rho"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#####################################################################
print(opt$file)
  name<-gsub("_dscaled.tab", "", opt$file)
  print(paste("name:",name))
  
  rho<-opt$rho
  print(paste("rho:",rho))
  
  M<-read.table(opt$file, sep="\t", header=T); 
  print(paste("dim M:",dim(M)))
  
  Cdouble <- CRobCor(M)
  gListDouble<-CreateGraph(Cdouble,rho) 
  gList<-gListDouble
  g<-gList[[1]]

  NetworkDesc(g, name,	Cdouble)


