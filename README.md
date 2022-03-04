# SpliceNet

To  systematically  explore  the  functions  of  core  splicing  factors  and regulators  in  alternative  splicing,  RNA-seq  analyses  were  carried  out  upon  the knockdown  of  over  304  splicing-related  factors.  

# Main scripts

Scripts associated to transcriptome-wide analysis of the effects of systematic knock down of splicing factors and regulators using siRNAs in HeLa cells.
All scripts were written in R. 

## Network FDR computation

Requiered libraries:
* R version >4.0
* parallel
* Hmisc
* ibble
* tidyr
* utils
* dplyr

## Prepare tables
PSI calculation was performed using VAST-TOOLS (jjjj9) version...
RAWDTA:

## Scaling

## Network FDR
fdr_CL_cor.R (use by CL or interactively)

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

Usage using CL: Rscript Net-fdr_CL_cor.R -f sscaled.tab  -s 0.1 -e 0.2 -i 0.05  -r 5 -n short -b ../SpliceNet 

* Output:
The function returns a table and their corresponding plots for the number of links for real and random data in the interval of correlations
FDR is computed as Number of Links RANDOM data / Number of links in REAL data * 100


For our files (see files):

* Net-fdr_CL_cor.R -s 0.1 -e 0.8 -i 0.01 -r 1000 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet -c 12     -f  A3_all_sscaled.tab -n A3long
* Net-fdr_CL_cor.R -s 0.1 -e 0.8 -i 0.01 -r 1000 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet -c 12  -f A5_all_sscaled.tab -n A5long
* Net-fdr_CL_cor.R -s 0.1 -e 0.8 -i 0.01 -r 1000 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet -c 12  -f ES_all_sscaled.tab -n ESlong
* Net-fdr_CL_cor.R -s 0.1 -e 0.8 -i 0.01 -r 1000 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet -c 12  -f IR_all_sscaled.tab -n IRlong

## Single Cor
* singlecor.R (use by CL or interactively)

Compute single correlation for a given dPSI table
* Input: Numeric matrix (only dPSI values, scaled or not), mininum correlation value, scripts folder, sample name

- Usage:

* Output: for a given threshold of correlation, it returns the edgelist: 
-   so, tg: source and target original order
-   Pearson cor and absCor, with  pvalue and fdr,
-  source and target alphabetically ordered,
-   link's name

* For our files:

* singlecor.R -m 0 -f A3_all_sscaled.tab  -p  1 -n A3 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
* singlecor.R -m 0 -f A5_all_sscaled.tab  -p 1 -n A5 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
* singlecor.R -m 0 -f ES_all_sscaled.tab  -p 1 -n ES -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
* singlecor.R -m 0 -f IR_all_sscaled.tab  -p 1 -n IR -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet

# Extract centrality
-  Net-centrality_CL_cor.R (for command line)
-  Net-centrality_interactive_cor.R (to run interactively)

For a given input matrix, it extracts degree (and normalized degree) poe each KD in a given interval of correlations
It can help to identify high/low stable factors.

#  Replace NAs: Analyze NAs in VAST-TOOLS INCLUSION final Table: 

* We count how many NAs contain a quality score for every PSI value computed. If they have more than 3, we replace this PSI value with NAnew3. 
* Input is the default vast-tools table (INCLUSION-...)
* Output: desired output, script 1: only replacing NAnew3

* If INCLUSION table was produced using last version of vast-tools, we can replace Intron Retention PSI value with NAI  if the adjusted pvalue for the "read umbalance" is below 0.05. 

* Input is the default vast-tools table (INCLUSION-...)
* Output: desired output, script 1: replacing NAnew3+NAI

## How to run:
- 1.- Download scripts in your desired folder.
- 2.- Load libraries:

library(data.table)
library(stringr)

- 3.- Source the function you want to use:

- source("replaceNA.R") # only NewNAs
- source("replaceNAI.R") # NewNAs + NAIs

* Define the input / output files:

- INCFile<-"test.tab"
- OUTFile<-"test_NA_NewN3.tab"
- OUTFile2<-"test_NA_NewN3_NAIs.tab"

- replaceNA(INCFile, OUTFile)
- replaceNAI(INCFile, OUTFile2)

