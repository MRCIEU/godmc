#!/bin/bash

set -e
source config
exec &> >(tee ${section_01_logfile})
print_version

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
	${phenotype_list} \
	${quality_scores} \
	${quality_scores_plot} \
	${methylation_summary}

# Check missingness
${plink} --bfile ${bfile_raw} --missing gz --out ${section_01_dir}/data

nrow=`zcat ./results/02/data.imiss.gz |awk 'NR>1 && $4==0 {print $0}'  |wc -l`
if [ ! "${nrow}" != "0" ]
then
	echo "Error: Your genotype data contains missing values. Please don't use a genotype probability cut-off."
    exit 1
else
    echo "You successfully converted genotype data to bestguess format without any missing values"
fi


