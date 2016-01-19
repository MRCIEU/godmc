#!/bin/bash

set -e
source config
exec &> >(tee ${section_10_logfile})
print_version


${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${smoking_pred}.plink \
	--qcovar ${gwas_covariates}.smoking \
	--out ${section_10_dir}/smoking \
	--thread-num ${nthreads}

echo "Compressing results"
gzip -f ${section_10_dir}/smoking.loco.mlma


echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
	${section_10_dir}/smoking.loco.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	0 \
	0 \
	0 \
	0 \
	${section_10_dir}/smoking.loco.mlma


echo "Successfully performed GWAS"
