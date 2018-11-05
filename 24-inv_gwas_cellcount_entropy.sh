#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_24_logfile})
print_version

if [ ! -f ${cellcounts_gwa} ]
then
    	echo "No multivariate cellcounts GWA will be performed; please note we only run GWAS on Houseman estimates"
        exit 0
fi

age=`awk '{print $4}' <${gwas_covariates}.cellcounts.numeric |sort -u |wc -l`
echo "Age variable has $age levels"

sex=`awk '{print $3}' <${gwas_covariates}.cellcounts.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${section_24_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--autosome \
	--out ${section_24_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_tf}.entropy.plink \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${section_24_dir}/cellcount_entropy \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_24_dir}/cellcount_entropy.loco.mlma


echo "Successfully performed GWAS"
