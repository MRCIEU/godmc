#!/bin/bash

set -e
source config
exec &> >(tee ${snp_data_logfile})


# Provision for having dosage data and convert to best guess if the analyst doesn't have this

# Copy raw data to modifiable data


# rewrite this to use keep ids in common with meth
echo "Copying genetic data to processing folder"
# cp ${bfile_raw}.bed ${bfile}.bed
# cp ${bfile_raw}.bim ${bfile}.bim
# cp ${bfile_raw}.fam ${bfile}.fam

${plink} \
	--bfile ${bfile_raw} \
	--keep ${intersect_ids_plink} \
	--maf ${snp_maf} \
	--hwe ${snp_hwe} \
	--geno ${snp_miss} \
	--mind ${snp_imiss} \
	--make-bed \
	--autosome \
	--out ${bfile} \
	--threads ${nthreads}


# Change SNP ids to chr:position:{SNP/INDEL}
echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (length($5) == "1" && length($6) == "1") print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${bfile}.bim

# Checking for any duplicate SNPs
cp ${bfile}.bim ${bfile}.bim.original2
awk '{
	if (++dup[$2] > 1) { 
		print $1, $2".duplicate."dup[$2], $3, $4, $5, $6 
	} else {
		print $0 
	}}' ${bfile}.bim.original2 > ${bfile}.bim

grep "duplicate" ${bfile}.bim | awk '{ print $2 }' > ${bfile}.duplicates.txt
plink1.90 \
	--bfile ${bfile} \
	--exclude ${bfile}.duplicates.txt \
	--make-bed \
	--out ${bfile} \
	--threads ${nthreads}

rm ${bfile}.*~

# Get frequencies
plink1.90 --bfile ${bfile} --freq --hardy --out ${matrixeqtl_mqtl_dir}/data
gzip -f ${matrixeqtl_mqtl_dir}/data.frq
gzip -f ${matrixeqtl_mqtl_dir}/data.hwe
gzip -c ${quality_scores} > ${matrixeqtl_mqtl_dir}/data.info.gz


# Make GRMs
echo "Creating kinship matrix"
gunzip -c ${hm3_snps} > temp_hm3snps.txt
${plink} \
	--bfile ${bfile} \
	--extract temp_hm3snps.txt \
	--maf ${grm_maf_cutoff} \
	--make-grm-bin \
	--out ${grmfile_all} \
	--threads ${nthreads}
rm temp_hm3snps.txt

# Create pedigree matrix if family data, otherwise remove related individuals from existing kinship and data file
if [ "${unrelated}" = "no" ]
then
	echo "Creating pedigree GRM"
	Rscript resources/relateds/grm_relateds.R ${grmfile_all} ${grmfile_relateds} ${rel_cutoff}
elif [ "${unrelated}" = "yes" ]
then
	echo "Removing any cryptic relateds"
	${plink} \
		--grm-bin ${grmfile_all} \
		--rel-cutoff ${rel_cutoff} \
		--make-grm-bin \
		--out ${grmfile_unrelateds} \
		--threads ${nthreads}

	${plink}  \
		--bfile ${bfile} \
		--keep ${grmfile_unrelateds}.grm.id \
		--make-bed \
		--out ${bfile} \
		--threads ${nthreads}
else 
	echo "Error: Set unrelated flag in config to yes or no"
	exit 1
fi

#Calculate PCs
gunzip -c ${hm3_snps_no_ld} > temp_hm3snpsnold.txt

${plink} \
	--bfile ${bfile} \
	--extract temp_hm3snpsnold.txt \
	--indep-pairwise 10000 5 0.1 \
	--maf 0.2 \
	--out ${pca} \
	--threads ${nthreads}

rm temp_hm3snpsnold.txt

if [ "${unrelated}" = "yes" ]
then
	${plink} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--pca 20 \
		--out ${pca} \
		--threads ${nthreads}
else

	${plink} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--make-bed \
		--out ${bfile}_ldpruned \
		--threads ${nthreads}

	Rscript resources/genetics/pcs_relateds.R \
		${bfile}_ldpruned \
		${pca} \
		${n_pcs} \
		${nthreads}
fi


# Get genetic outliers
echo "Detecting genetic outliers"

Rscript resources/genetics/genetic_outliers.R \
	${pcs_all} \
	${pca_sd} \
	${n_pcs} \
	${genetic_outlier_ids} \
	${pcaplot}

# Remove genetic outliers from data
echo "Removing genetic outliers"
${plink} \
	--bfile ${bfile} \
	--remove ${genetic_outlier_ids} \
	--make-bed \
	--out ${bfile} \
	--threads ${nthreads}

${gcta} \
	--grm ${grmfile_all} \
	--remove ${genetic_outlier_ids} \
	--make-grm-bin \
	--out ${grmfile_all} \
	--thread-num ${nthreads}



#Update ids
awk '{print $1,$2}' <${bfile}.fam >${intersect_ids_plink}
awk '{print $2}' <${bfile}.fam >${intersect_ids}

echo "Successfully formatted SNP data"