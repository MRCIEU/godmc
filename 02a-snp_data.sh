#!/bin/bash

set -e
source config
exec &> >(tee ${section_02a_logfile})
print_version

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
${plink} \
	--bfile ${bfile} \
	--exclude ${bfile}.duplicates.txt \
	--make-bed \
	--out ${bfile} \
	--threads ${nthreads}


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
if [ "${related}" = "yes" ]
then
	echo "Creating pedigree GRM"
	Rscript resources/relateds/grm_relateds.R ${grmfile_all} ${grmfile_relateds} ${rel_cutoff}
elif [ "${related}" = "no" ]
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
	echo "Error: Set related flag in config to yes or no"
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

if [ "${related}" = "no" ]
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



# If there are any genetic outliers then remove them and recalculate PCs
# Otherwise don't do anything

n_outliers=`wc -l ${genetic_outlier_ids} | awk '{ print $1 }'`
if [ "${n_outliers}" = "0" ]
then
	echo "No genetic outliers detected"
else 
	# Remove genetic outliers from data
	echo "Removing ${n_outliers} genetic outliers from data"
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
#From here on, we have clean data
	echo "Recalculating PCs with outliers removed"

	if [ "${related}" = "no" ]
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

fi

#Recode alleles
cp ${bfile}.bim ${bfile}.bim.original3

Rscript resources/genetics/harmonization.R \
	${bfile}.bim \
	${SNPfail_allelecoding}
	 

# Get frequencies
${plink} --bfile ${bfile} --freq gz --hardy gz --missing gz --out ${section_02_dir}/data

gzip -f -c ${quality_scores} > ${section_02_dir}/data.info.gz

nrow=`zcat ./results/02/data.imiss.gz |awk 'NR>1 && $4==0 {print $0}'  |wc -l`
if [ ! "${nrow}" != "0" ]
then
	echo "Error: Your genotype data contains missing values. Please don't use a genotype probability cut-off."
    exit 1
else
    echo "You successfully converted genotype data to bestguess format without any missing values"
fi


#Update ids
awk '{print $1,$2}' < ${bfile}.fam > ${intersect_ids_plink}
awk '{print $2}' < ${bfile}.fam > ${intersect_ids}

rm -f ${bfile}.*~
rm temp_hm3snpsnold.txt


echo "Successfully formatted SNP data"