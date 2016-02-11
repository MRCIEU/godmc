#!/bin/bash

#Height prediction
set -e
source ./config
exec &> >(tee ${section_03b_logfile})
print_version

if [ ! "${phenotypes}" = "NULL" ]
then
	echo "Calculating height allele score"

	${plink} \
		--bfile ${bfile} \
		--keep ${intersect_ids_plink} \
		--extract ${height_snps_wood}.snp \
		--indep-pairwise 10000 5 0.1 \
		--out ${height_snps}

	fgrep -v -f ${height_snps}.prune.out ${height_snps_wood}.txt > ${height_snps_wood}.pruned

	${plink} \
		--bfile ${bfile} \
		--keep ${intersect_ids_plink} \
		--score ${height_snps_wood}.pruned \
		--out ${height_snps}

	# var EA beta
	Rscript ./resources/genetics/calculateheightcor.R \
		${height_snps}.profile \
		${ewastransformed} \
		${height_plot}
	
	echo "Successfully completed allele score calculation "
else
	echo "No phenotypes"
fi