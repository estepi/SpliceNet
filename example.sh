
#!/usr/bin/bash


module use /software/as/el7.2/EasyBuild/CRG/modules/all 
module load R/4.1.0-foss-2021a 

export R_LIBS=/users/jvalcarcel/emancini/R/packages-4.1.0

cd /no_backup/jvalcarcel/SpliceNet


R CMD Rscript prepareALLTable.R 
R CMD Rscript prepareExonsTable.R 
R CMD Rscript prepareIRTable.R 
R CMD Rscript prepareA5Table.R 
R CMD Rscript prepareA3Table.R 


R CMD Rscript prepareExonsTable_subset3500.R 
R CMD Rscript prepareIRTable_subset3500.R 
R CMD Rscript prepareA5Table_subset3500.R 
R CMD Rscript prepareA3Table_subset3500.R 



