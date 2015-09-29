#!/bin/bash

set -e
source config


# Provision for having dosage data and convert to best guess if the analyst doesn't have this

# Copy raw data to modifiable data


# rewrite this to use keep ids in common with meth
echo "Copying genetic data to processing folder"
cp ${bfile_raw}.bed ${bfile}.bed
cp ${bfile_raw}.bim ${bfile}.bim
cp ${bfile_raw}.fam ${bfile}.fam

# Change SNP ids to chr:position:{SNP/INDEL}
echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (length($5) == "1" && length($6) == "1") print "chr"$1, $1":"$4":SNP", $3, $4, $5, $6;else print "chr"$1, $1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${bfile}.bim


# Make GRMs and calculate PCs
echo "Creating kinship matrix and calculating PCs"
gunzip -c ${hm3_snps} > temp_hm3snps.txt
${plink} \
	--bfile ${bfile} \
	--extract temp_hm3snps.txt \
	--maf ${grm_maf_cutoff} \
	--pca ${n_pcs} \
	--make-grm-bin \
	--out ${grmfile_all}
rm temp_hm3snps.txt

# Get genetic outliers
echo "Detecting genetic outliers"
R --no-save --args ${pcs_all} ${pca_sd} ${n_pcs} ${genetic_outlier_ids} < resources/genetics/genetic_outliers.R

# Remove genetic outliers from data
echo "Removing genetic outliers"
${plink} \
	--bfile ${bfile} \
	--remove ${genetic_outlier_ids} \
	--make-bed \
	--out ${bfile}

${gcta} \
	--grm ${grmfile_all} \
	--remove ${genetic_outlier_ids} \
	--make-grm-bin \
	--out ${grmfile_all}

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
		--out ${grmfile_unrelateds}

	${plink}  \
		--bfile ${bfile} \
		--keep ${grmfile_unrelateds}.grm.id \
		--make-bed \
		--out ${bfile}
else 
	echo "Error: Set unrelated flag in config to yes or no"
	exit 1
fi
