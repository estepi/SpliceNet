# SpliceNet

To  systematically  explore  the  functions  of  core  splicing  factors  and regulators  in  alternative  splicing,  RNA-seq  analyses  were  carried  out  upon  the knockdown  of  over  304  splicing-related  factors.  

# Main scripts

Scripts associated to transcriptome-wide analysis of the effects of systematic knock down of splicing factors and regulators using siRNAs in HeLa cells.
All scripts were written in R. 

## Network FDR computation

Requiered libraries:
* R version: 
* optparse
* parallel
* tibble
* tidyr
* utils


## Prepare tables
PSI calculation was performed using VAST-TOOLS (jjjj9) version...
RAWDTA:

## Scaling

## fdr_CL_cor.R
This script computes Network FDR represented as the ratio betwenn TRUE links and RANDOM links. Input file is a matrix with EVENTS in rows and KDs (samples) in the column

dPSI values can be none, single or double scaled (should be prepared in advance, see :...)

The function retunrs the total number of links which are below a given correlation value and are significant (corrected pvalue <0.1).
Expected values are higher FDR at lower correlation values

Matrix randomization is prepared randomizing rows, so correlation between columns are disrupted
As noiser is the data, more random links will be find, higher FDR

* Parameters:

- -s 0.1 (start: from which correlation value the function scan the data. Correlation values are between 0 and 1)
- -e 0.4  (end: till which correlation value the funcion scan the data. Correlation values are between 0 and 1)
- -i 0.02 (interval: size of the interval to compute the correlation.  For example, 0.02 means you will scan 0, 0.02, 0.04, 0.06, etc)
- -r 100 (number of random matrixes. Nuber of links in random data come from the mean of the number of links of those matrixs)
- -b (bin: folder where functions_fdr_MC_cor.R is)
- -c 12 (cores: number of cores to use)
- -f sscaled.tab (file: input file, should be prepared in advance. It contains only dPSI values, scaled or not  ex: sscaled.tab) 
- -n A3short (name: prefix to use for output files, ej: A3short)

* Output:
The function returns a table and their corresponding plots for the number of links for real and random data in the interval of correlations
FDR is computed as Number of Links RANDOM data / Number of links in REAL data * 100

## Single Cor
* singlecor.R

Compute single correlation for a given dPSI table
* Input: Numeric matrix (only dPSI values, scaled or not), mininum correlation value, scripts folder, sample name
 
* Output: for a given threshold of correlation, it returns the edgelist: 
+   so, tg: source and target original order
+   Pearson cor and absCor, with  pvalue and fdr,
+  source and target alphabetically ordered,
+   link's name

# Extract centrality
-  Net-centrality_CL_cor.R (for command line)
-  Net-centrality_interactive_cor.R (to run interactively)

For a given input matrix, it extracts degree (and normalized degree) poe each KD in a given interval of correlations
It can help to identify high/low stable factors.
