library(optparse)
library(utils)
library(igraph)
library(scales)
source("functions_fdr_noMC_cor.R")
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


}

#####################################################################
option_list = list(
   make_option(c("-f", "--file"), type="character", default=NULL, 
              help="eventscaled file name", metavar="character"),
    make_option(c("-m", "--min"), type="numeric", default="0.1",
              help="min cor value"),
    make_option(c("-f", "--fdr"), type="numeric", default="1",
              help="max fdr value"),
    make_option(c("-n", "--name"), type="character", default="test"))
  

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

M<-read.table(opt$file, sep="\t", header=T); 
print(paste("dim M:",dim(M)))

MCor <- CRobCor(M)
g<-graph_from_data_frame(MCor[,1:2])
NetworkDesc(g, name,	MCor)


