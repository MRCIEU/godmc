#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_22_logfile})
print_version

#age=`awk '{print $3}' <${gwas_covariates}.aar.numeric |sort -u |wc -l`
#echo "Age variable has $age levels"
if [ ! -f ${gwas_covariates}.aar.numeric ]
then
	echo "No Age acceleration GWA will be performed as age is not variable"
	exit 0
fi

sex=`awk '{print $3}' <${gwas_covariates}.aar.factor |sort -u |wc -l`
if [ "$sex" -gt "0" ]
then

${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar.numeric \
	--covar ${gwas_covariates}.aar.factor \
	--autosome \
	--out ${section_22_dir}/aar \
	--thread-num ${nthreads}
fi


if [ "$sex" -eq "1" ]
then

${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${age_pred}.aar.plink \
	--qcovar ${gwas_covariates}.aar.numeric \
	--autosome \
	--out ${section_22_dir}/aar \
	--thread-num ${nthreads}
fi


echo "Compressing results"
gzip -f ${section_22_dir}/aar.loco.mlma

echo "Successfully performed IWAS"
