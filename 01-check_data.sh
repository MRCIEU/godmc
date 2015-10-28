#!/bin/bash

set -e
source config

# R --no-save --args ${bfile_raw}.bim ${bfile_raw}.fam ${betas} ${covariates} ${provided_cellcounts} ${intersect_ids} ${intersect_ids_plink} ${snpcheck} ${snpplot} ${poscheck} ${EWASphenotypes} ${height} ${bmi} ${cnvs} ${log_directory} < resources/datacheck/datacheck.R

Rscript resources/datacheck/datacheck.R ${bfile_raw}.bim ${bfile_raw}.fam ${betas} ${covariates} ${provided_cellcounts} ${intersect_ids} ${intersect_ids_plink} ${snpcheck} ${control_snps} ${phenotypes} ${cnvs} ${cohort_descriptives}
