#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_22f_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else
	echo ""
	echo "Performing BMI IWAS"
	echo ""

  Rscript resources/inversions/adapt_pheno.R ${ewastransformed} BMI ${inv_processed_dir}/

	${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${inv_processed_dir}/BMI.plink \
	--autosome \
	--out ${section_22_dir}/BMI \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_22_dir}/BMI.loco.mlma

echo "Successfully performed BMI EWAS"
