###################################################
#functions
flat_cor_mat <- function(cor_r, cor_p) {
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor,-1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p,-1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}
###################################################
set.seed(1234)
createRandomMatrix <- function(TrueM, NumRandomM)
{
  if (is.numeric(TrueM))
  {
    MList <- vector("list", NumRandomM)
    for (i in 1:NumRandomM) {
#    MList[[i]] <- t(apply(TrueM, 1, sample, replace=TRUE))
     MList[[i]] <- t(apply(TrueM, 1, sample, replace=FALSE))
     colnames(MList[[i]])<-colnames(TrueM)

    }
    
      return(MList)  }
  else
  {
    print("You are not entering a Numerical matrix")
  }
  
}
####################################################################
getEdges <- function(C, rho, fdr)
{
  
  df <-    rcorr(C)
  cor_r <- df$r
  cor_p <- df$P
  my_cor_matrix <- flat_cor_mat(cor_r, cor_p)
  my_cor_matrixClean <-
  my_cor_matrix[my_cor_matrix$row != my_cor_matrix$column,]
  my_cor_matrixClean$fdr <- p.adjust(my_cor_matrixClean$p, method = "fdr")
  links <- my_cor_matrixClean[my_cor_matrixClean$cor > rho &
                             my_cor_matrixClean$fdr < fdr, ]
  print(length( E(g) ))
  print(dim(links))
  print(head(links))
  return(nrow(links))
}

getCorM <- function(C, rho, fdr)
{
  
  df <-    rcorr(Ms)
  cor_r <- df$r
  cor_p <- df$P
  
  df2 <- flat_cor_mat(cor_r, cor_p)
  colnames(df2)[1:2]<-c("so","tg")
  
  df2Clean <-
    df2[df2$so != df2$tg,]
    df2Clean$absCor<-abs( df2Clean$cor)  
    df2Clean$fdr <- p.adjust(df2Clean$p, method = "fdr")
  
  #order links alphabetically and put rownames:
  edgeListLCSort<-data.frame(t(unlist(apply(df2Clean[,1:2],
                                            1, function(x){sort(x)}))))
  final<-data.frame(df2Clean,"source"=  edgeListLCSort$X1, 
                    "target" =edgeListLCSort$X2,
                    link=paste(edgeListLCSort$X1,
                               edgeListLCSort$X2, sep="-")
                    )
  

  return(final)
}

######################################################################

#agregar from to como parametros
getEdgesByRHO <- function(RList, start, end, interval, ncores, fdr)
{
  MList <- vector("list", length(seq(start, end, interval)))
  for (rho in seq(start,
                  end, interval))
  {
    res <- unlist(mclapply(RList,  getEdges,  rho,  mc.cores = ncores, fdr=fdr))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}
######################################################################
#rho sera el valor de corte de la correlacion ahora
getEdgesBySample <- function(rrlist, rho, fdr)
{
  
  RRList <- vector("list")
  i = 1
for (i in 1:length(rrlist))
  {
    res <- mean(unlist(lapply(rrlist[[i]],  getEdges,  rho, fdr)))
    RRList[[i]] <- res
    i <- i + 1
  }
  return(RRList)
}
#############################################################
#esta es la funcion que vamos a modificar poniendo el threshold en la correlacion

getEdgesByRHO <- function(RList, start, end, interval, ncores, fdr)
{
  MList <- vector("list", length(seq(start, end, interval)))
  i = 1
  
  for (rho in seq(start,
                  end, interval))
  {
    res <- unlist(mclapply(RList,  
                           getEdges, 
                           rho,  mc.cores = ncores, fdr))
    MList[[i]] <- res
    i <- i + 1
  }
  return(MList)
}
#########################Plots###########################################
plots <- function(start,
                  end,
                  interval,
                  trueEdges,
                  randomEdges,
                 name)
{
  ratioNumberOfEdges <- randomEdges / trueEdges * 100
  eje1 <- seq(start, end, interval)
  name <- name
  ###########################################################################################
  file1 <- paste(paste(paste(
    "NumberOfEdges_", paste(start, end, sep = "-"), sep = ""
  )
  , name, sep = "_"),
  ".png", sep = "")
  png(file1)
  plot(
    eje1,
    trueEdges,
    main = "Number of edges real / random data vs rho" ,
    type = "l",
    xlab = "rho",
    ylab = "Number of edges",
    col = "blue",
    ylim = c(min(trueEdges, randomEdges), max(trueEdges, randomEdges)),
    xlim = c(start, end)
  )
  lines(eje1, randomEdges, col = "red")
  legend(
    "topright",
    c("real data", "random data"),
    lty = 1,
    col = c("blue", "red")
  )
  dev.off()
  message("Plot 1: Number of edges real / random data vs rho")
  ###########################################################################################
  file1df <- paste(paste(paste(
    "NumberOfEdges_", paste(start, end, sep = "-"), sep = ""
  )
  , name, sep = "_"),
  ".tab", sep = "")
  df1 <-
    data.frame(
      rho = eje1,
      random = randomEdges,
      true = trueEdges,
      ratio = randomEdges / trueEdges
    )
  write.table(df1,
              file = file1df,
              sep = "\t",
              col.names = NA)
  ###########################################################################################
  file2 <- paste(paste(paste(
    "RatioNumberOfEdges_",
    paste(start, end, sep = "-"), sep = ""
  ),
  name, sep = "_"),
  ".png", sep = "")
  png(file2)
  plot(eje1,
       (randomEdges / trueEdges) * 100,
       main = "Ratio edges in random data vs real data",
       type = "l", main=name)
  abline(h = c(0, 5, 10))
  dev.off()
  message("Plot 2: Ratio edges in random data vs real data")
}
###########################################################
estimateFDR<-function(
  inputM=sampleData,#eventscaled
  NumRandomM=10,
  start=0.3,
  end=1,
  interval=0.05,
  ncores=1,
  name,
  fdr=0.05)
{
  TrueList<-vector("list",1)
  TrueList[[1]]<-sampleData
  TPList<-rep(TrueList, NumRandomM)
  TrueEdgesList<-getEdgesByRHO(TPList, start, end, interval, ncores, fdr)
  
  print(paste("Lenght of true edge list:", length(TrueEdgesList)   )  ) 
  RRList<-createRandomMatrix(TrueM = sampleData,  NumRandomM)
  RandomEdgesList<-getEdgesByRHO(RRList, start, end, interval, ncores, fdr)
  
  print(paste("Lenght of random edge list:",length(RandomEdgesList))) 
  trueEdges<-unlist(lapply(TrueEdgesList, mean))
  randomEdges<- unlist(lapply(RandomEdgesList, mean))
##################################################
  #add an small correcction in order to avoid infinite values
  message("Finish computation. Lets plot.")
  
print( start)
print( end)
print(interval)
print(trueEdges+0.01)
print(randomEdges+0.01)
print(name)#all together

plots( start,  end, interval, 
        trueEdges=trueEdges+0.01,
        randomEdges=randomEdges+0.01,
        name)#all together
##################################################
}



  
