#!/bin/bash

set -e
source config
exec &> >(tee ${section_01_logfile})
print_version

Rscript resources/datacheck/requirements.R

Rscript resources/datacheck/genetic_data.R \
	${bfile_raw}.bim \
	${bfile_raw}.fam \
	${quality_scores} \
	${control_snps} \
	${snpchrtxt} \
	${snpchrplot} \
	${genetic_descriptives} \
	${quality_scores_plot}

Rscript resources/datacheck/methylation_data.R \
	${betas} \
	${bfile_raw}.fam \
	${provided_cellcounts} \
	${meth_ids} \
	${methylation_descriptives} \
	${methylation_summary} \
	${intersect_ids} \
	${intersect_ids_plink}

Rscript resources/datacheck/covariates.R \
	${covariates} \
	${bfile_raw}.fam \
	${meth_ids} \
	${ageplot} \
	${covariate_descriptives}

Rscript resources/datacheck/phenotypes.R \
	${phenotypes} \
	${meth_ids} \
	${covariates} \
	${phenotype_descriptives} \
	${phenotype_list}

Rscript resources/datacheck/cnv_data.R \
	${cnvs} \
	${meth_ids} \
	${covariates} \
	${cnv_descriptives}

Rscript resources/datacheck/collect_descriptives.R \
	${genetic_descriptives} \
	${methylation_descriptives} \
	${covariate_descriptives} \
	${phenotype_descriptives} \
	${cnv_descriptives}	\
	${cohort_descriptives}

echo "You successfully performed all data checks"


