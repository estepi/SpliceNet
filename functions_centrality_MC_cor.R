###################################################
#functions
flat_cor_mat <- function(cor_r, cor_p) {
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor,-1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p,-1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}
########################################################
getCorM <- function(M, min, fdr)
{
  df <-    rcorr(M)
  cor_r <- df$r
  cor_p <- df$P
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
#OK,extract netwotks from the different MATRIXS
######################################################################

getCentralityByCOR <- function(M, start, end, interval, fdr)
{
  vrange<-seq(start,  end, interval)
  resMatrixDG<-data.frame(matrix(ncol=length(vrange), nrow = ncol(M)))
  colnames(resMatrixDG)<-vrange
  rownames(resMatrixDG)<-colnames(M)
  head(resMatrixDG)
  
  resMatrixNDG<-data.frame(matrix(ncol=length(vrange), nrow = ncol(M)))
  colnames(resMatrixNDG)<-vrange
  rownames(resMatrixNDG)<-colnames(M)
  head(resMatrixNDG)
  
  MList <- vector("list", length(vrange))
  DGList <- vector("list", length(vrange))
  NDGList <- vector("list", length(vrange))
  
  
  i=1
  for (i in 1:length(vrange))
    {
    print(i)
    print(vrange[i])
    res <- getCorM(sampleData, vrange[i],  fdr=fdr)
    g<-graph_from_data_frame(res[, 7:8], directed = FALSE)
    DG <- degree(g)
    NDG <- DG / length(E(g))#
    
    MList[[i]] <- res
    DGList[[i]] <- DG
    NDGList[[i]] <- NDG
    i <- i + 1
  }
  
  p=1
  for(p in 1:length(vrange))
  {
    dd <- match(rownames(resMatrixDG), names(unlist(DGList[[p]])))
    resMatrixDG[, p] <- unlist(DGList[[p]])[dd]
    
    nn <- match(rownames(resMatrixNDG), names(unlist(NDGList[[p]])))
    resMatrixNDG[, p] <- unlist(NDGList[[p]])[nn]
  }
  
  ResList<-list(resMatrixDG, resMatrixNDG)
  return(ResList)
}


