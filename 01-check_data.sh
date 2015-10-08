#!/bin/bash

set -e
source config

R --no-save --args ${bfile_raw}.bim ${bfile_raw}.fam ${betas} ${covariates} ${provided_cellcounts} ${home_directory}/processed_data/ids/ids.txt ${home_directory}/processed_data/ids/ids_plink.txt ${home_directory}/processed_data/datacheck/no.snps.by.chr ${home_directory}/processed_data/datacheck/no.snps.by.chr.plot ${home_directory}/resources/genetics/snpsforbuildandposcheck.txt ${EWASphenotypes} ${height} ${bmi} < resources/datacheck/datacheck.R ${cnvs} > ${home_directory}/processed_data/datacheck/datacheck.Rout 2>&1

which git
which PLINK
which GCTA
