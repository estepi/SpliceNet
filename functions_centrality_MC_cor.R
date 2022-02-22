###################################################
#functions
flat_cor_mat <- function(cor_r, cor_p) {
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor,-1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p,-1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}
####################################################################
getEdges <- function(M, min, fdr)
  # return edges
{
  df <-    rcorr(M)
  cor_r <- df$r
  cor_p <- df$P
  summary(cor_p)
  df2 <- flat_cor_mat(cor_r, cor_p)
  colnames(df2)[1:2]<-c("so","tg")
  head(df2)
  df2Clean <-   df2[df2$so != df2$tg,]
  df2Clean$absCor<-abs( df2Clean$cor)  
  df2Clean$p
  df2Clean$fdr <- p.adjust(df2Clean$p, method = "fdr")
  head(df2Clean$fdr)
  
  links <- df2Clean[df2Clean$absCor > min &
                      df2Clean$fdr < fdr,]
  
  edgeListLCSort<-data.frame(t(unlist(apply(links[,1:2],
                                            1, function(x){sort(x)}))))
  head(edgeListLCSort)
  final<-data.frame(links,"source"=  edgeListLCSort$X1, 
                    "target" =edgeListLCSort$X2,
                     link=paste(edgeListLCSort$X1,
                               edgeListLCSort$X2, sep="-"))
  finalUnique<-final[!duplicated(final$link),]
  print(head(finalUnique))
  return(nrow(finalUnique))
}
########################################################
getCorM <- function(M, min, fdr)
{
  df <-    rcorr(M)
  head(df)
  
  cor_r <- df$r
  head(cor_r)
  
  cor_p <- df$P
  head(cor_p)
  
  df2 <- flat_cor_mat(cor_r, cor_p)
  colnames(df2)[1:2]<-c("so","tg")
  
  
  df2Clean <-
  df2[df2$so != df2$tg,]
  df2Clean$absCor<-abs( df2Clean$cor)  
  df2Clean$fdr <- p.adjust(df2Clean$p, method = "fdr")
  
  #filter here
  links <- df2Clean[df2Clean$absCor > min &
                      df2Clean$fdr < fdr,]
  
  #order links alphabetically and put rownames:
  edgeListLCSort<-data.frame(t(unlist(apply(links[,1:2],
                                        1, function(x){sort(x)}))))
  
  final<-data.frame(links,"source"=  edgeListLCSort$X1, 
                    "target" =edgeListLCSort$X2,
                    link=paste(edgeListLCSort$X1,
                               edgeListLCSort$X2, sep="-"))
  #remove duplicates
  finalUnique<-final[!duplicated(final$link),]
  print(head(finalUnique))
  return(finalUnique)
  
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
computeCentrality<-function(
  inputM=sampleData,#eventscaled
  start=0.3,
  end=1,
  interval=0.05,
  ncores=1)
{
  TrueList<-vector("list",1)
  TrueList[[1]]<-sampleData
#cargao la matrix   
  #input matriz output conectivity dif correlation values
TrueEdgesList<-getEdgesByRHO(TPList, start, end, interval, ncores)

message("Finish computation")
  

}



NetworkDesc <- function(g)
{
###################################################
  DG <- degree(g)
  NDG <- DG / length(E(g))#
}
