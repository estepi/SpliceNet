subsetFromGCDistribution <- function(dPSIFile = table,
                                     all = exons_gc,
                                     name,
                                     outputPath) {
  ###############################################################
  #    dPSIFile <- table
  #   head(dPSIFile)
  #    all <- exons_gc
  #    name
  #    outputPath
  ##############################################
  file1 <- paste(outputPath,
                 paste(name, "GCDist.png", sep = "."), sep = "/")
  png(
    file1,
    width = 12,
    height = 8,
    units = 'in',
    res = 300
  )
  hist(
    dPSIFile$GC,
    breaks = 100,
    xlim = c(0, 1) ,
    main = paste(name, "GC distribution", sep = " ")
  )
  dev.off()
  
  ##############################################
  
  file2 <- paste(outputPath,
                 paste(name, "allExons_GCDist.png", sep = "."),
                 sep = "/")
  png(
    file2,
    width = 12,
    height = 8,
    units = 'in',
    res = 300
  )
  hist(
    all$EXON_GCC,
    breaks = 100,
    xlim = c(0, 1) ,
    main = paste(name, "All exons GC distribution", sep = " ")
  )
  dev.off()
  
  ##############################################
  file3 <-
    paste(outputPath, paste(name, "boxplots.png", sep = "."), sep = "/")
  png(
    file3,
    width = 12,
    height = 8,
    units = 'in',
    res = 300
  )
  boxplot(dPSIFile$GC, main = "GC distribution")
  abline(h = c(boxplot(dPSIFile$GC)$stats))
  dev.off()
  ##############################################
  file4 <-
    paste(outputPath,
          paste(name, "allExons_boxplots.png", sep = "."),
          sep = "/")
  png(
    file4,
    width = 12,
    height = 8,
    units = 'in',
    res = 300
  )
  boxplot(all$EXON_GCC, main = "All exons GC distribution")
  abline(h = c(boxplot(all$EXON_GCC)$stats))
  dev.off()
  
  ##############################################
  #recomputo las distribuciones:
  Q1L <- boxplot(dPSIFile$GC)$stats[2,]
  MedianL <- boxplot(dPSIFile$GC)$stats[3,]
  Q3L <- boxplot(dPSIFile$GC)$stats[4,]
  ######################################################################
  ######################################################################
  print(paste("GC Q1 value ", round(Q1L, digits = 2), sep = ": "))
  print(paste("Num of Events lower than GC Q1 ",
              length(which(dPSIFile$GC <= Q1L)), sep = ": "))
  
  
  #single scaled
  Q1dPSI <- dPSIFile[dPSIFile$GC <= Q1L,]
  Q1dPSI$AVcontrols <- NULL
  Q1dPSI$GC <- NULL
  dim(Q1dPSI)#4169*305
  
  ##################################################
  outQ1_GC <-
    paste(outputPath, paste(name, "Q1_GC_dPSI.tab", sep = "_"), sep = "/")
  write.table(Q1dPSI, outQ1_GC,  sep = "\t")
  ##################################################
  Q1Numeric <- data.frame(matrix(
    unlist(apply(Q1dPSI, 1, as.numeric)),
    ncol = ncol(Q1dPSI),
    nrow = nrow(Q1dPSI),
    byrow = T
  ))
  colnames(Q1Numeric) <- colnames(Q1dPSI)
  rownames(Q1Numeric) <- rownames(Q1dPSI)
  
  deltascaled <- scale(Q1Numeric)# scaled by columns (KDs)
  
  outQ1_GCss <-
    paste(outputPath,
          paste(name, "Q1_GC_dPSI_sscaled.tab", sep = "_"),
          sep = "/")
  
  write.table(deltascaled, outQ1_GCss,  sep = "\t")
  #double scaled
  dsscaled <- t(scale(t(deltascaled)))    #scaled by events
  
  outQ1_GCds <-
    paste(outputPath,
          paste(name, "Q1_GC_dPSI_dscaled.tab", sep = "_"),
          sep = "/")
  
  write.table(dsscaled, outQ1_GCds,  sep = "\t")
  
  
  Q1Events <- rownames(dPSIFile)[dPSIFile$GC <= Q1L]
  ######################################################################
  print(paste("GC Q3 value ", Q3L, sep = ": "))
  print(paste("Num of Events greater than GC Q3 ",
              length(which(dPSIFile$GC >= Q3L)), sep = ": "))
  
  
  #single scaled
  Q3dPSI <- dPSIFile[dPSIFile$GC >= Q3L,]
  Q3dPSI$AVcontrols <- NULL
  Q3dPSI$GC <- NULL
  dim(Q3dPSI)#4168*305
  
  ##################################################
  outQ3_GC <-
    
    paste(outputPath, paste(name, "Q3_GC_dPSI.tab", sep = "_"), sep = "/")
  write.table(Q3dPSI, outQ1_GC,  sep = "\t")
  ##################################################
  Q3Numeric <- data.frame(matrix(
    unlist(apply(Q3dPSI, 1, as.numeric)),
    ncol = ncol(Q3dPSI),
    nrow = nrow(Q3dPSI),
    byrow = T
  ))
  colnames(Q3Numeric) <- colnames(Q3dPSI)
  rownames(Q3Numeric) <- rownames(Q3dPSI)
  
  deltascaled <- scale(Q3Numeric)# scaled by columns (KDs)
  
  outQ3_GCss <-
    paste(outputPath,
          paste(name, "Q3_GC_dPSI_sscaled.tab", sep = "_"),
          sep = "/")
  
  write.table(deltascaled, outQ1_GCss,  sep = "\t")
  #double scaled
  dsscaled <- t(scale(t(deltascaled)))    #scaled by events
  
  outQ3_GCds <-
    paste(outputPath,
          paste(name, "Q3_GC_dPSI_dscaled.tab", sep = "_"),
          sep = "/")
  
  write.table(dsscaled, outQ3_GCds,  sep = "\t")
  
  Q3Events <- rownames(dPSIFile)[dPSIFile$GC >= Q3L]
  
  ###############################################
  #events are classify in HIGH, LOwer and MID GC
  
  ee <- match(rownames(dPSIFile), rownames(all))
  ddDF <- all[ee, ]
  aggDFEvents <- data.frame(aggregate(
    ddDF$EVENT,
    by = list(ddDF$GENEID),
    FUN = paste,
    collapse = ":"
  ))
  colnames(aggDFEvents) <- c("gene", "events")
  head(aggDFEvents)
  
  dim(aggDFEvents)#4115 genes
  Q1i <- match(Q1Events, rownames(ddDF))
  Q3i <- match(Q3Events, rownames(ddDF))
  
  ddDF$GC_content <- rep("MED", nrow(ddDF))
  ddDF$GC_content[Q1i] <- "LOW"
  ddDF$GC_content[Q3i] <- "HIGH"
  table(ddDF$GC_content)
  
  
  #sort by gene and PRINT:
  out_GC <-
    paste(outputPath, paste(name, "gene_GC.tab", sep = "_"), sep = "/")
  ddDF_Sorted <- ddDF[order(ddDF$GENEID), ]
  head(ddDF_Sorted)
  write.table(ddDF_Sorted,
              file = out_GC,
              col.names = NA,
              sep = "\t")
  #Agrregate by GENE
  aggDF <-
    data.frame(aggregate(
      ddDF_Sorted$GC_content,
      by = list(ddDF$GENEID),
      FUN = paste,
      collapse = ":"
    ))
  colnames(aggDF) <- c("gene", "gcContent")
  head(aggDF)
  
  out_GCAgg <-
    paste(outputPath, paste(name, "gene_Agg_GC.tab", sep = "_"), sep = "/")
  write.table(aggDF,
              file = out_GCAgg,
              sep = "\t",
              col.names = NA)
  print(paste("Files printed in", outputPath, sep = ": "))
  
  
}
###############################################################################
