# Transcriptome-wide analysis of the effects of systematic knock down of splicing factors and regulators using siRNAs in HeLa cells

## Description
RNA-seq upon knockdown of 304 splicing-related factors and regulators using siRNAs in HeLa cells

### Sample collection: 
HeLa cells were collected after 72hr post transfection. 

### Treatment protocol:
HeLa cells were forward transfected in 96-well plates with a siRNA library against known factors and regulators using SMARTpools (ON TARGET PLUS, Thermo-Scientific) using an automatized robotic procedure (Sciclone Liquid Handling Workstation, Perkin Elmer). 

### Growth protocol 
HeLa were purchased from the American Type Culture Collection (ATCC). Cells were cultivated in Glutamax Dulbecco’s modified Eagle’s medium (Life technologies) supplemented with 10% Fetal bovine serum (Life technologies) and antibiotics (penicillin 500 u/ml; streptomycin 0.5 mg/ml, Life Technologies).

### Nucleic Acid extraction protocol 
Endogenous	mRNAs were purified 72 hours post	 transfection by using	 oligo	 dT-coated 96 well plates	(mRNA catcher PLUS, Life Technologies) following the manufacturer’s instructions.	 

### Nucleic Acid library construction protocol
Ultra low input library preparation. Full-length RNA-seq from single cells using Smart-seq2

### Nucleic Acid Sequencing protocol
Illumina HiSeq 2000 / Illumina HiSeq 2500 2x125 bp

* Sequencing quality check was done using FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
* Raw Data is available here: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11202

### HTS Alignment protocol

Read alignment and gene count were performed using STAR version 2.5.3a gainst Homo sapiens genome assembly v38 (ENSEMBL) 
Gene count tables were prepared using STAR option --quantMode GeneCounts. It includes only uniquely mapped reads in the quantification.
* Pipelines automatization: NEXTFLOW  ~  version 0.27.1
* Read alignmenNot default parameters:
                --outSJfilterReads Unique 
                --outFilterType BySJout 
                --outFilterMultimapNmax 10 
                --alignSJoverhangMin 6 
                --alignSJDBoverhangMin 3 
                --outFilterMismatchNoverLmax 0.1 
                --alignIntronMin 20 
        	      --alignIntronMax 500 
	              --outSAMstrandField intronMotif 
                --outFilterIntronMotifs RemoveNoncanonicalUnannotated 
                --seedSearchStartLmax 50 

* Alignments were quality checked with QUALIMAP. 

### CLASS COLORS class_colors_2020.txt
* Factors were classified into  Family Class according XX
* (Assembbly) Order was defined according


### GE profiles: Gene_logFC_GW*.tab

* Differential gene expression analysis was performed using egdeR R/Bioconductor package (http://bioconductor.org/packages/release/bioc/html/edgeR.html).

* Low expressed genes were filtered: at least 70% of the samples should have more than 1 cpm. We include knock down’s gene, even if its expression is lower than 1 cpm. Total genes analyzed: **11876**

### AS profiles

* Alternative splicing genome wide quantification was performed using vast-tools (https://github.com/vastgroup/vast-tools) against Homo sapiens version 38 (ensembl) (HS2).**CHECK VERSION**

* Quality checked was performed and those  samples/genes/ with… Nas were discarded **(Check with GOSIA)**
We add an extra layer of filtering discarding thos IR events with *unbalanced* coverage (pval <0.1)

* Alterantive Splicing event type which are analyzed are: ES, IR, A5, A3

* We consider an event differentially spliced if:
- there were not Nas in the row (?)
- standard deviation of controls was lower than 5 (?)
- delta PSI was greater than 25 (?)

### Total initial input (dPSI_full_No_Nas.txt): 
* 82795 events well mapped
* Filter by RANGE > 5, total: 41870

 - **Alt3**: 6093

- **Alt5**: 3676

- **Exons**: 16721
  
        *  **ANN**: 1088
	*  **C1**: 1513
	*  **C2**: 1359
	*  **C3**: 1772
	*  **MIC**: 72
	*  **S**: 10917 
	
 - **IR**:  15380

## Prepare dPSI, sscaled and dscaled files for main types of EVENTS: ES, A5, A3 and IR

- Scripts are: **prepareEVENTtable.R**

Using /SpliceNetData//class_colors_2020.txt and XXX under SpliceNetData/ folder, yo can get 3 files from each event type:

* EVENT_dPSI.tab
* EVENT_all_sscaled.tab (data scaled by columns)
* EVENT_all_dscaled.tab (data double scaled)



## AS Cross Regulatory Network: Corss_edges.csv
 From AS profile table, we only consider events which belong to Knock Down genes and  dpSI > 25. **CHECK GOSIA** 
 
 From this analysis we can estimate the effect of cross regulation between splicing factors via alternative splicing
 
 In this table you can filter by Source, Target and dPSI
 
 Total edges: **9966 **
 
 ### Total input:


### AS Functional Network

In order to decipher functional relationships between splicing factors, we computed Pearson pairwise correlations
For each dataset (event type, ES, IR, A5, A3) we computed Network FDR and  Pearson correlation and corresponding pvalue for each pair. We report all possible correlations

 * Total expected links: 46360 (*(305*305)-305*))

* Input:

     - A3_all_sscaled.tab: 6094 events
     - A5_all_sscaled.tab: : 3677 events
     - ES_all_sscaled.tab: 16722 events
     - IR_all_sscaled.tab: 15381 events 
     
* Script parameter	

- singlecor.R -m 0 -f A3_all_sscaled.tab  -p  1 -n A3 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
- singlecor.R -m 0 -f A5_all_sscaled.tab  -p 1 -n A5 -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
- singlecor.R -m 0 -f ES_all_sscaled.tab  -p 1 -n ES -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet
- singlecor.R -m 0 -f IR_all_sscaled.tab  -p 1 -n IR -b /no_backup/jvalcarcel/emancini/SpliceNet/SpliceNet

### 
use values form F scores (GOSIA will wrtie and exlain)
0.28 for general + 

