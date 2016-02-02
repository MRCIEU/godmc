#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_11_logfile})
print_version

age=`awk '{print $4}' <${gwas_covariates}.cellcounts.numeric |sort -u |wc -l`
echo "Age variable has $age levels"

sex=`awk '{print $3}' <${gwas_covariates}.cellcounts.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--covar ${gwas_covariates}.cellcounts.factor \
	--out ${section_11_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--out ${section_11_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--covar ${gwas_covariates}.cellcounts.factor \
	--out ${section_11_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_11_dir}/cellcount_entropy.loco.mlma


echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
	${section_11_dir}/cellcount_entropy.loco.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	0 \
	0 \
	0 \
	0 \
	${section_11_dir}/cellcount_entropy.loco.mlma


echo "Successfully performed GWAS"
