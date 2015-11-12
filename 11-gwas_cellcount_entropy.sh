#!/bin/bash

set -e
source config
exec &> >(tee ${section_11_logfile})
print_version


${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts \
	--out ${section_11_dir}/cellcount_entropy \
	--thread-num ${nthreads}

echo "Compressing results"
gzip -f ${section_11_dir}/cellcount_entropy.loco.mlma

echo "Successfully performed GWAS"
