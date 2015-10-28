#!/bin/bash

set -e
source config

Rscript resources/datacheck/datacheck.R \
	${bfile_raw}.bim \
	${bfile_raw}.fam \
	${betas} \
	${covariates} \
	${provided_cellcounts} \
	${intersect_ids} \
	${intersect_ids_plink} \
	${snpcheck} \
	${snpplot} \
	${control_snps} \
	${phenotypes} \
	${cnvs} \
	${cohort_descriptives} \
	${ageplot} \
	2>&1 | tee ${datacheck_logfile}
