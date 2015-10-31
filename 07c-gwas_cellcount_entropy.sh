#!/bin/bash

set -e
source config
exec &> >(tee ${gwas_cellcount_entropy_logfile})

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts \
	--out ${gwas_cellcounts_dir}/cellcount_entropy \
	--thread-num ${nthreads}

echo "Compressing results"
gzip -f ${gwas_cellcounts_dir}/cellcount_entropy.loco.mlma

echo "Successfully performed GWAS"
