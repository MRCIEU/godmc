#!/bin/bash

set -e
source config


# Provision for having dosage data and convert to best guess if the analyst doesn't have this

# Copy raw data to modifiable data

echo "Copying genetic data to processing folder"
cp ${bfile_raw}.bed ${bfile}.bed
cp ${bfile_raw}.bim ${bfile}.bim
cp ${bfile_raw}.fam ${bfile}.fam

# Change SNP ids to chr:position:{SNP/INDEL}
echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (length($5) == "1" && length($6) == "1") print $1, $1":"$4":SNP", $3, $4, $5, $6;else print $1, $1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${bfile}.bim


# Make GRMs and calculate PCs
echo "Creating kinship matrix and calculating PCs"
gunzip -c ${hm3_snps} > temp_hm3snps.txt
${plink} \
	--bfile ${bfile} \
	--extract temp_hm3snps.txt \
	--maf ${grm_maf_cutoff} \
	--pca ${npcs} \
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








# Generate meQTL format

#Generate plink.raw / plink.frq #freq files are required to determine effect allele
echo "Converting plink files to transposed raw format"
${plink} --bfile ${bfile} --recode A-transpose --out ${bfile} --freq

echo "Generating matrixeqtl format genetic data"
cut -f 2,5,6 ${bfile}.traw > ${allele_ref}

head -n 1 ${bfile}.traw | cut -f 7- | tr '\t' '\n' | tr '_' '\t' | cut -f 2 | tr '\n' '\t' | awk '{ printf ("snpid\t%s\n", $0) }' > traw.header

sed 1d ${bfile}.traw | cut -f 2,7- > ${bfile}.traw2
mv ${bfile}.traw2 ${bfile}.traw

echo "Splitting data into ${bfile_chunksize} SNP subsets"
split -d -a 10 -l ${bfile_chunksize} ${bfile}.traw ${tabfile}.tab.

i=1
for file in ${tabfile}.tab.*
do
	echo ${file}
	cat traw.header ${file} > ${tabfile}.tab.${i}
	rm ${file}
	i=$(($i+1))
done

rm traw.header


echo "Done!"

#covariates file #ID COV1 COV2

#matrixQTL
#id indiv1 indiv1
#gender 0 1
#age 36 41

#genotype
#id indiv1 indiv2
#snp1 2 0
#snp2 0 1

#id indiv1 indiv2
#CpG1 0.21 0.11
#CpG2 0.1 0.16

#geneloc chr s1 s2
#CpG1 chr1 1 100
#CpG2 chr1 200 300

#snp loc
#snp chr pos
#snp1 chr1 1
#snp2 chr1 5



# ./resources/genetics/make_matrixeqtl_format.py ${bfile} ${bfile_chunksize}

# This will make the files
# ${bfile}.1.tab, ${bfile}.2.tab, ${bfile}.3.tab etc

# make GRM from hm3 SNPs
# if family data then make pedigree matrix
# if unrelateds then calculate principal components
# convert bim file to change rs IDs to chr<x>:<pos>
# remove cryptically related individuals #still to do

# create matrix eqtl format 
## hash has this already use --recode option

#To do:
#what about PCs on related data? Do you correct for pop strat by using kinship matrix?
#what about covariates?
#remove related indiv
#check SNPs on right build
#check allele coding, no R, ID, D anything other than ACTG
