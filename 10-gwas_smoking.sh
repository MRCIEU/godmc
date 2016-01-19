#!/bin/bash

set -e
source ./config

declare -a pheno=('' ${smoking_pred} ${smoking_pred}_ge25 ${smoking_pred}_lt25)

batch=${1}

if [[ "$batch" -lt "1" ]] || [[ "$batch" -gt "3" ]]
then
	echo "Error - please specify a batch variable between 1 and 3. e.g."
	echo "./10-gwas_smoking.sh 1"
	exit
fi

pheno_file=${pheno[${batch}]}
pheno=`basename ${pheno_file}`

exec &> >(tee ${section_10_logfile})
print_version


${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${pheno_file}.plink \
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
