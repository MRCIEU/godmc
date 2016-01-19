#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_03a_logfile})
print_version

if [ ! "${phenotypes}" = "NULL" ]
then
	echo "Transforming phenotypes"
	Rscript	resources/phenotypes/phenotype_transform.R \
		${intersect_ids_plink} \
		${phenotypes} \
		${ewasplot} \
		${covariates} \
		${ewastransformed} \
		${ewas_sd} \
		${ewas_plink} \
		${ewas_file}
	echo "Successfully completed transformation"
else
	echo "No phenotypes"
fi

