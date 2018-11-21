#!/bin/bash

set -e
source ./config

# declare -a pheno=('' ${smoking_pred} ${smoking_pred}_ge25 ${smoking_pred}_lt25)

IFS=$'\r\n' GLOBIGNORE='*' command eval  'pheno=($(cat ${section_10_dir}/gwas_list.txt))'
ngwas=${#pheno[@]}
batch=${1}
if [[ "$batch" -lt "1" ]] || [[ "$batch" -gt "${ngwas}" ]]
then
	echo "Error - please specify a batch variable between 1 and ${ngwas}"
	echo "Usage: ${0} [batch]"
	exit
fi

batch=$(( ${1} - 1 ))

echo "The following GWASs can be performed:"
for i in $(seq 1 $ngwas)
do
	ii=$(( $i - 1 ))
	echo "${i}. ${pheno[${ii}]}"
done

echo ""
echo "Performing:"
echo "${1}. ${pheno[${batch}]}"

pheno_file="${methylation_processed_dir}/${pheno[${batch}]}.plink"
pheno=${pheno[${batch}]}
cov=`echo ${pheno/#smoking_prediction/}`
exec &> >(tee ${section_23_logfile}${1})
print_version

age=`awk '{print $3}' <processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric |sort -u |wc -l`
echo "Age variable has $age levels"
sex=`awk '{print $3}' <processed_data/genetic_data/gwas_covariates${cov}.smoking.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--autosome \
	--out ${section_23_dir}/${pheno} \
	--thread-num ${nthreads}
fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--autosome \
	--out ${section_23_dir}/${pheno} \
	--thread-num ${nthreads}
fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then

${gcta} \
	--bfile ${bfile_inv} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--autosome \
	--out ${section_23_dir}/${pheno} \
	--thread-num ${nthreads}
fi

echo "Compressing results"
gzip -f ${section_23_dir}/${pheno}.loco.mlma

echo "Successfully performed GWAS"