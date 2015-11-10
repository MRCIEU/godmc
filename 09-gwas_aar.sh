#!/bin/bash

set -e
source config
exec &> >(tee ${section_09_logfile})

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar \
	--out ${section_09_dir}/aar \
	--thread-num ${nthreads}

echo "Compressing results"
gzip -f ${section_09_dir}/aar.loco.mlma

echo "Successfully performed GWAS"