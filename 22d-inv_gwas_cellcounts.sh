#!/bin/bash

set -e
source ./config

if [ ! -f ${cellcounts_gwa} ]
then
    	echo "No multivariate cellcounts GWA will be performed; please note we only run GWAS on Houseman estimates"
        exit 0
fi

batch=${1}
re='^[0-9]+$'
ncellcounts=`wc -l ${section_22_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
if ! [[ $batch =~ $re ]] ; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

if [ "${batch}" -gt "${ncellcounts}" ]; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

if [ "${batch}" -lt "1" ]; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

exec &> >(tee ${section_22d_logfile}${batch})
print_version

age=`awk '{print $4}' <${gwas_covariates}.cellcounts.numeric |sort -u |wc -l`
echo "Age variable has $age levels"

sex=`awk '{print $3}' <${gwas_covariates}.cellcounts.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_plink} \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${section_22_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_plink} \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--autosome \
	--out ${section_22_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${cellcounts_plink} \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${section_22_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi


echo "Compressing results"
gzip -f ${section_22_dir}/cellcount_${batch}.loco.mlma


echo "Successfully performed GWAS for cell type ${batch}"
