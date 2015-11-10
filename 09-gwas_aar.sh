#!/bin/bash

set -e
source config
exec &> >(tee ${section_09_logfile})

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar \
	--out ${gwas_aar_dir}/aar \
	--thread-num ${nthreads}

echo "Compressing results"
gzip -f ${gwas_aar_dir}/aar.loco.mlma

echo "Successfully performed GWAS"