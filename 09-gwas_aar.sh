#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_09_logfile})
print_version

age=`awk '{print $3}' <${gwas_covariates}.aar.numeric |sort -u |wc -l`
echo "Age variable has $age levels"
sex=`awk '{print $3}' <${gwas_covariates}.aar.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar.numeric \
	--covar ${gwas_covariates}.aar.factor \
	--out ${section_09_dir}/aar \
	--thread-num ${nthreads}
fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar.numeric \
	--out ${section_09_dir}/aar \
	--thread-num ${nthreads}
fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--covar ${gwas_covariates}.aar.factor \
	--out ${section_09_dir}/aar \
	--thread-num ${nthreads}
fi


echo "Compressing results"
gzip -f ${section_09_dir}/aar.loco.mlma


echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
	${section_09_dir}/aar.loco.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	0 \
	0 \
	0 \
	0 \
	${section_09_dir}/aar.loco.mlma


echo "Successfully performed GWAS"