#!/bin/bash

set -e
source ./config
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
	--out ${bfile} \
	--chr 1-23 \
	--threads ${nthreads}


# Sex check

n23=`grep ^23 ${bfile}.bim | wc -l`

if [ "$n23" -gt "0" ]
then
	
	${plink} \
		--bfile ${bfile} \
		--split-x b37 no-fail \
		--make-bed \
		--out ${bfile}

	${plink} \
		--bfile ${bfile} \
		--check-sex \
		--out ${section_02_dir}/data
	
	nprob=`grep "PROBLEM" ${section_02_dir}/data.sexcheck |wc -l`

    if [ "$nprob" -eq "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
	fi


	if [ "$nprob" -gt "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
		echo "They will be removed."
		echo "The summary is located here:"
		echo "${section_02_dir}/data.sexcheck"

		grep "PROBLEM" ${section_02_dir}/data.sexcheck | awk '{print $1, $2}' > ${bfile}.failed_sexcheck

		${plink} \
			--bfile ${bfile} \
			--remove ${bfile}.failed_sexcheck \
			--make-bed \
			--out ${bfile}
	fi

fi

# Change SNP ids to chr:position:{SNP/INDEL}
echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") &&  ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G")) print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${bfile}.bim

#Recode alleles to uniform format eg. I/D for INDELs
cp ${bfile}.bim ${bfile}.bim.original2
touch ${SNPfail1}
touch ${bfile}.duplicates.txt

Rscript resources/genetics/harmonization.R \
	${bfile}.bim \
	${SNPfail1}

# Checking for any duplicate SNPs
cp ${bfile}.bim ${bfile}.bim.original3
awk '{
	if (++dup[$2] > 1) { 
		print $1, $2".duplicate."dup[$2], $3, $4, $5, $6 
	} else {
		print $0 
	}}' ${bfile}.bim.original3 > ${bfile}.bim

grep "duplicate" ${bfile}.bim | awk '{ print $2 }' > ${bfile}.duplicates.txt

#Remove SNPs with low info scores
awk '$3 < 0.80 {print $1}' <${quality_scores} > ${bfile}.lowinfoSNPs.txt

cat ${bfile}.duplicates.txt ${SNPfail_allelecoding} ${bfile}.lowinfoSNPs.txt |sort -u >${bfile}.failed.SNPs.txt

n_failedSNPs=`wc -l ${bfile}.failed.SNPs.txt | awk '{ print $1 }'`

# Remove SNPs from data
echo "Removing ${n_failedSNPs} SNPs from data"


${plink} \
	--bfile ${bfile} \
	--exclude ${bfile}.failed.SNPs.txt \
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
	--threads ${nthreads} \
	--autosome
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
	--autosome \
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

fi


# Find mismatched SNPs and misaligned SNPs with EasyQC
# Get frequencis for strand check
${plink} \
	--bfile ${bfile} \
	--freq \
	--out ${bfile}

if [ -f "processed_data/genetic_data/easyQC.ecf.out" ]
then
	echo "easyqc files present from previous run which will be removed"
	rm processed_data/genetic_data/easy*
	rm processed_data/genetic_data/*easy*
else
	echo "passed file check"
fi

    
Rscript ./resources/genetics/easyQC.R ${bfile}.bim ${bfile}.frq ${easyQC} ${easyQCfile} ${easyQCscript}
mv ./processed_data/genetic_data/easyQC.multi.AFCHECK.png ./results/02
mv ./processed_data/genetic_data/easyQC.rep ./results/02

# Remove mismatched SNPs and flip misaligned SNPs
echo "Remove mismatched SNPs and NO FLIPPING"

${plink} \
	--bfile ${bfile} \
	--exclude ${easyQC}.mismatch_afcheck.failed.SNPs.txt \
	--make-bed \
	--out ${bfile} \
	--threads ${nthreads}


# From here on, we have clean data

if [ ! "${n_outliers}" = "0" ]
then

	echo "Recalculating PCs with outliers removed"

	if [ "${related}" = "no" ]
	then
		${plink} \
			--bfile ${bfile} \
			--extract ${pca}.prune.in \
			--pca 20 \
			--out ${pca} \
			--autosome \
			--threads ${nthreads}
	else

		${plink} \
			--bfile ${bfile} \
			--extract ${pca}.prune.in \
			--make-bed \
			--out ${bfile}_ldpruned \
			--autosome \
			--threads ${nthreads}

		Rscript resources/genetics/pcs_relateds.R \
			${bfile}_ldpruned \
			${pca} \
			${n_pcs} \
			${nthreads}
	fi

fi

# Get frequencies, missingness, hwe, info scores
${plink} \
	--bfile ${bfile} \
	--freq gz \
	--hardy gz \
	--missing gz \
	--out ${section_02_dir}/data

gzip -f -c ${quality_scores} > ${section_02_dir}/data.info.gz

# Check missingness
missingness=`zcat ${section_02_dir}/data.imiss | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'`

echo "Average missingness: ${missingness}"

if (( $(bc <<< "${missingness} > 0.02") ))
then
	echo ""
	echo ""
	echo ""
	echo ""
	echo "WARNING"
	echo ""
	echo ""
	echo "Your genetic data has missingness of ${missingness}"
	echo ""
	echo "This seems high considering that you should have converted to best guess format with a very high hard call threshold"
	echo ""
	echo "Please ensure that this has been done"
fi

# Update ids
awk '{print $1,$2}' < ${bfile}.fam > ${intersect_ids_plink}
awk '{print $2}' < ${bfile}.fam > ${intersect_ids}

rm -f ${bfile}.*~
rm temp_hm3snpsnold.txt


echo "Successfully formatted SNP data"
