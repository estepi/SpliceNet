# SpliceNet
Scripts associated to transcriptome-wide analysis of the effects of systematic knock down of splicing factors and regulators using siRNAs in HeLa cells

# Description

To  systematically  explore  the  functions  of  core  splicing  factors  and regulators  in  alternative  splicing,  RNA-seq  analyses  were  carried  out  upon  the knockdown  of  over  304  splicing-related  factors.  

# Network FDR computation
Install libraries:
* optparse
* parallel
* tibble
* tidyr
* utils


# Prepare tables

# Single Cor

## fdr_CL_cor.R
This script computes Network FDR represented as the ratio betwenn TRUE links and RANDOM links. Input file is a matrix with EVENTS in rows and KDs (samples) in the column

dPSI values can be none, single or double scaled (should be prepared in advance, see :...)

The function retunrs the total number of links which are below a given correlation value and are significant (pvalue <0.1).
Expected values are higher FDR at lower correlation values

Matrix randomization is prepared randomizing rows, so correlation between columns are disrupted
As noiser is the data, more random links will be find, higher FDR

## Parameters:

- -s 0.1 (start: from which correlation value the function scan the data. Correlation values are between 0 and 1)
- -e 0.4  (end: till which correlation value the funcion scan the data. Correlation values are between 0 and 1)
- -i 0.02 (interval: size of the interval to compute the correlation.  For example, 0.02 means you will scan 0, 0.02, 0.04, 0.06, etc)
- -r 100 (number of random matrixes. Nuber of links in random data come from the mean of the number of links of those matrixs)
- -b (bin: folder where the scripts are)
- -c 12 (cores: number of cores to use)
- -f sscaled.tab (file: input file, should be prepared in advance. It contains only dPSI values, scaled or not  ex: sscaled.tab) 
- -n A3short (name: prefix to use for output files, ej: A3short)

# Extract centrality

