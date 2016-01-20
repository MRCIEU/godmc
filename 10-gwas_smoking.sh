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

exec &> >(tee ${section_10_logfile}${batch})
print_version


${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${pheno_file} \
	--qcovar ${gwas_covariates}.smoking.numeric \
	--covar ${gwas_covariates}.smoking.factor \
	--out ${section_10_dir}/${pheno} \
	--thread-num ${nthreads}

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
