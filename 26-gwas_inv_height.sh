#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_26_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else
	echo ""
	echo "Performing Height IWAS"
	echo ""

  Rscript resources/inversions/adapt_pheno.R ${ewastransformed} ${ewas_plink} Height ${inv_processed_dir}/

	${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${inv_processed_dir}/Height.plink \
	--autosome \
	--out ${section_26_dir}/Height \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_26_dir}/Height.loco.mlma

echo "Successfully performed Height EWAS"
