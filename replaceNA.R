replaceNA <- function(INCFile,
                      OUTFile) {
  #######################################################################
  #READ INPUT FILE
  INCData <- read.csv(
    INCFile,
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    row.names = 2
  )
  
  #######################################################################
  #subset PSI values
  onlyValues <- INCData[, 6:ncol(INCData)]
  #subset QUAL values
  Qvals <- grep("\\.Q$", colnames(onlyValues))
  qualValues <- onlyValues[, Qvals]
  psiValues  <- onlyValues[, -Qvals]
  finalValues <- matrix(0, ncol = ncol(psiValues), nrow = nrow(psiValues))
  colnames(finalValues) <- colnames(onlyValues)[-Qvals]
  #######################################################################
  #count Ns in QUAL values
  Ns_comma <-
    t(apply(qualValues, 1, function(x) {
      str_count(x, "N,")
    }))
  #######################################################################
  #replace
  indexes_N <- which(Ns_comma != 0, arr.in = TRUE)
  indexes_NewN3 <- which(Ns_comma == 3, arr.in = TRUE)
  
  finalValues[indexes_N] <- Ns_comma[indexes_N]
  head(finalValues)
  correctedPsiVal <- psiValues
  
  
  correctedPsiVal[indexes_NewN3] <- "NAnew3"
  originalNA <- is.na(psiValues)
  correctedPsiVal[originalNA] <- "NAold"
  fullNAS <- cbind(INCData[, 1:5], correctedPsiVal)
  #####################################################
  #compute and print STATS:
  tVal <- length(originalNA)# 235328422
  oNA <- length(which(originalNA))#95995460 40.8
  new3 <- length(indexes_NewN3) / 2# 827226 28.39
  oNA2 <- length(which(correctedPsiVal == "NAold"))#
  new3b <- length(which(correctedPsiVal == "NAnew3"))#
  
  m1 <- paste(round(oNA / tVal * 100), "%")
  m2 <- paste(round(new3 / tVal * 100), "%")
  m3 <- paste(round(new3b / tVal * 100), "%")
  
  message(paste("Num total PSI values:", tVal))
  message(paste("\n NA original values:",  m1))
  message(paste("\n NewN3 values:", m3))
  message(paste("\n Ttoal NA + NewN3 values:", m2))
  #####################################################
  #PRINT OUTPUT FILE
  write.table(fullNAS, OUTFile, sep = "\t", col.names = NA)
  message(paste("\n Final corrected table is written in:", OUTFile))
}