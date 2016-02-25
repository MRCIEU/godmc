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

pheno_file="${methylation_processed_dir}/${pheno[${batch}]}.plink"
pheno=${pheno[${batch}]}
cov=`echo ${pheno/#smoking_prediction/}`
exec &> >(tee ${section_10_logfile}${1})
print_version

age=`awk '{print $3}' <processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric |sort -u |wc -l`
echo "Age variable has $age levels"
sex=`awk '{print $3}' <processed_data/genetic_data/gwas_covariates${cov}.smoking.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

n23=`grep ^23 ${bfile}.bim | wc -l`

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--autosome \
	--out ${section_10_dir}/${pheno} \
	--thread-num ${nthreads}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--grm ${grmfile_all} \
	--out ${section_10_dir}/${pheno}_chr23 \
	--thread-num ${nthreads}

fi

if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--autosome \
	--out ${section_10_dir}/${pheno} \
	--thread-num ${nthreads}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${pheno_file} \
	--qcovar processed_data/genetic_data/gwas_covariates${cov}.smoking.numeric \
	--grm ${grmfile_all} \
	--out ${section_10_dir}/${pheno}_chr23 \
	--thread-num ${nthreads}

fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--autosome \
	--out ${section_10_dir}/${pheno} \
	--thread-num ${nthreads}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${pheno_file} \
	--covar processed_data/genetic_data/gwas_covariates${cov}.smoking.factor \
	--grm ${grmfile_all} \
	--out ${section_10_dir}/${pheno}_chr23 \
	--thread-num ${nthreads}

fi

head -n1 ${section_10_dir}/${pheno}.loco.mlma >${section_10_dir}/${pheno}.loco
tail -q -n +2 ${section_10_dir}/${pheno}.loco.mlma ${section_10_dir}/${pheno}_chr23.mlma >>${section_10_dir}/${pheno}.loco
mv ${section_10_dir}/${pheno}.loco ${section_10_dir}/${pheno}.loco.mlma
rm ${section_10_dir}/${pheno}_chr23.mlma

echo "Compressing results"
gzip -f ${section_10_dir}/${pheno}.loco.mlma


echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
	${section_10_dir}/${pheno}.loco.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	0 \
	0 \
	0 \
	0 \
	${section_10_dir}/${pheno}.loco.mlma


echo "Successfully performed GWAS"
