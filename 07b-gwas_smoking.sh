#!/bin/bash

set -e
source config
exec &> >(tee ${gwas_smoking_logfile})

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${smoking_pred}.plink \
	--qcovar ${gwas_covariates}.smoking \
	--out ${gwas_smoking_dir}/smoking \
	--thread-num ${nthreads}

echo "Compressing results"
gzip ${gwas_smoking_dir}/smoking.loco.mlma

echo "Successfully performed GWAS"
