#!/bin/bash

set -e
source config
batch=${1}
re='^[0-9]+$'
ncellcounts=`wc -l ${section_12_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
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

exec &> >(tee ${section_12_logfile}${batch})

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${cellcounts_plink} \
	--qcovar ${gwas_covariates}.cellcounts \
	--out ${section_12_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}

echo "Compressing results"
gzip -f ${section_12_dir}/cellcount_${batch}.loco.mlma

echo "Successfully performed GWAS for cell type ${batch}"
