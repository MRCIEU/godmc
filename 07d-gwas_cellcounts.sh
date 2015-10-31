#!/bin/bash

set -e
source config

cd ${gwas_cellcounts_dir}

chr=${1}
re='^[0-9]+$'
if ! [[ $chr =~ $re ]] ; then
	echo "error: Chromosome variable is not a number"
	echo "Usage: ${0} [chr]"
	exit 1
fi
exec &> >(tee ${gwas_cellcounts_logfile}_${chr})

echo "Formatting data for chromosome ${chr}"
plink1.90 \
	--bfile ${bfile} \
	--chr ${chr} \
	--make-bed \
	--out ${bfile}_${chr}

cut -d " " -f 1-5 ${bfile}_${chr}.fam | tr ' ' '\t' > ${bfile}_${chr}.fam.temp
paste -d "\t" ${bfile}_${chr}.fam.temp ${cellcounts_tf}.gemma > ${bfile}_${chr}.fam

nval=`awk '{ print NR }' ${gwas_cellcounts_dir}/cellcounts_columns.txt | tr '\n' ' '`
echo "Performing multivariate LMM on ${nval} cellcount phenotypes"


${gemma} \
	-bfile ${bfile}_${chr} \
	-c ${gwas_covariates}.cellcounts.gemma \
	-k ${grmfile_all}.gemma \
	-n ${nval} \
	-lmm 4 \
	-o cellcounts_mvlmm_${chr}

echo "Compressing results"
gzip -f output/cellcounts_mvlmm_${chr}.assoc.txt
mv output/cellcounts_mvlmm_${chr}.* ${gwas_cellcounts_dir}

rm ${bfile}_${chr}*

echo "Successfully performed multivariate LMM on cellcounts for chromsome ${chr}"
