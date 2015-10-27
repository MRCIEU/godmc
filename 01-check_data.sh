#!/bin/bash

set -e
source config

R --no-save --args ${bfile_raw}.bim ${bfile_raw}.fam ${betas} ${covariates} ${provided_cellcounts} ${intersect_ids} ${intersect_ids_plink} ${snpcheck} ${snpplot} ${poscheck} ${EWASphenotypes} ${height} ${bmi} ${cnvs} ${log_directory} ${ageplot} < resources/datacheck/datacheck.R > ${home_directory}/processed_data/datacheck/datacheck.Rout 2>&1

