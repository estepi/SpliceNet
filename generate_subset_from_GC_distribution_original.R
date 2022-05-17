#################################################################
subsetFromGCDistribution <- function(inputFile, name, outputPath) {
  ###############################################################
  onlyValues <- inputFile <- table
  print(paste("PSI values dim", dim(onlyValues), sep = ": "))
  file1 <-
    paste(outputPath, paste(name, "boxplots.png", sep = "."), sep = "/")
  png(
    file1,
    width = 12,
    height = 8,
    units = 'in',
    res = 300
  )
  boxplot(onlyValues$GC, main = "GC distribution")
  abline(h = c(boxplot(Q3_GC_Df$GC)$stats))
  dev.off()
##############################################
  #recomputo las distribuciones:
  LWL <- boxplot(onlyValues$GC)$stats[1,]
  Q1L <- boxplot(onlyValues$GC)$stats[2,]
  MedianL <- boxplot(onlyValues$GC)$stats[3,]
  Q3L <- boxplot(onlyValues$GC)$stats[4,]
  UWL <- boxplot(onlyValues$GC)$stats[5,]
  ######################################################################
  length(which(onlyValues$GC >= LWL))#todos
  length(which(onlyValues$GC <= LWL))#2
  length(which(onlyValues$GC <= Q1L))#
  length(which(onlyValues$GC <= Q1L & onlyValues$GC >= LWL))#
  length(which(onlyValues$GC >= Q3L))#
  length(which(onlyValues$GC <= UWL))#
  length(which(onlyValues$GC >= UWL))#
  #####################################################################
  #printing
  print(paste("GC LW value ", round(LWL, digits = 2), sep = ": "))
  print(paste("Num of Events lower than GC Lw ",
              length(which(onlyValues$GC <= LWL)), sep = ": "))
  outLW_GC <-
    paste(outputPath, paste(name, "LW_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC <= LWL, ], outLW_GC, col.names = NA, sep =
                "\t")
  ######################################################################
  print(paste("GC Q1 value ", round(Q1L, digits = 2), sep = ": "))
  print(paste("Num of Events lower than GC Q1 ", length(which(onlyValues$GC <=
                                                                Q1L)), sep = ": "))
  outQ1_GC <-
    paste(outputPath, paste(name, "Q1_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC <= Q1L, ], outQ1_GC, col.names = NA, sep =
                "\t")
  ######################################################################
  print(paste(round(LWL, digits = 2), "< GC < ", round(Q1L, digits = 2), "bp ", sep =
                " "))
  print(paste(
    "Num of Events greater than LW and lower than GC Q1 ",
    length(which(onlyValues$GC <= Q1L) &
             length(which(onlyValues$GC >= LWL)))
    ,
    sep = ": "
  ))
  outLW_Q1_GC <-
    paste(outputPath, paste(name, "LW_Q1_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC <= Q1L & onlyValues$GC >= LWL, ]
              ,
              outLW_Q1_GC,
              col.names = NA,
              sep = "\t")
  ######################################################################
  print(paste("GC Q3 value ", Q3L, sep = ": "))
  print(paste("Num of Events greater than GC Q3 ", length(which(onlyValues$GC >=
                                                                  Q3L)), sep = ": "))
  outQ3_GC <-
    paste(outputPath, paste(name, "Q3_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC >= Q3L, ], outQ3_GC, col.names = NA, sep =
                "\t")
  ######################################################################
  print(paste("GC UW value", round(UWL, digits = 2), sep = ": "))
  print(paste("Num of Events greater than Range UW",
              length(which(onlyValues$GC >= UWL)), sep = ": "))
  outWD_GC <-
    paste(outputPath, paste(name, "UW_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC >= UWL, ] , outWD_GC, col.names = NA, sep =
                "\t")
  ######################################################################
  print(paste(round(Q3L, digits = 2), "< GC < ", round(UWL, digits = 2), "bp ", sep =
                " "))
  print(paste(
    "Num of Events greater than Q3 and lower than GC UW ",
    length(which(onlyValues$GC >= Q3L &
                   onlyValues$GC <= UWL)),
    sep = ": "
  ))
  outQ3_UW_GC <-
    paste(outputPath, paste(name, "Q3_UW_GC_psi.tab", sep = "_"), sep = "/")
  write.table(onlyValues[onlyValues$GC >= Q3L & onlyValues$GC <= UWL, ],
              outQ3_UW_GC,
              col.names = NA,
              sep = "\t")
  ######################################################################
  print(paste("Files printed in", outputPath, sep = ": "))
}
###############################################################################
workingDir<-"inputGC/"
outputPath<-"inputGC/surfaces/"
inputFile<-"ES_psi.tab"
name="ES_R"
setwd(workingDir)
#############
table <-
  read.table(
    inputFile,
    sep = "\t",
    dec = "," ,
    header = T,
    row.names = 1
  )
colnames(table)
print(paste("Total input", nrow(table), sep = ": "))
#########################################
features <- read.table("exons_features_hs2.tab", sep = "\t", header = T)
rownames(features) <- features$EVENT
exons_gc <- features[, c(1, 19)]
write.table(features[, c(1, 19)],
            "exons_GC.tab",
            col.names = NA,
            sep = "\t")
#########################################
Q3_GCi <- match(rownames(table), rownames(features))
Q3_GC_Df <- data.frame(EVENT = row.names(features)[Q3_GCi],
                       GC = features$EXON_GCC[Q3_GCi])
table$GC <- features$EXON_GCC[Q3_GCi]
length(which(is.na(table$GC)))
table <- table[!is.na(table$GC), ]
###############################################################################
subsetFromGCDistribution(inputFile, name, outputPath)
###############################################################################
file2 <- paste(outputPath,
               paste(name, "GCDist.png", sep = "."), sep = "/")
png(
  file2,
  width = 12,
  height = 8,
  units = 'in',
  res = 300
)
hist(
  table$GC,
  breaks = 100,
  xlim = c(0, 1) ,
  main = paste(name, "GC distribution", sep = " ")
)
dev.off()
getwd()