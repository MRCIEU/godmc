#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_27_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else
	echo ""
	echo "Performing BMI IWAS"
	echo ""

  Rscript resources/inversions/adapt_pheno.R ${ewastransformed} ${ewas_plink} BMI ${inv_processed_dir}/

	${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${inv_processed_dir}/BMI.plink \
	--autosome \
	--out ${section_27_dir}/BMI \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_27_dir}/BMI.loco.mlma

echo "Successfully performed BMI EWAS"
