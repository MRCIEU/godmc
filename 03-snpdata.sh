#!/bin/bash

set -e
source config

# Provision for having dosage data and convert to best guess if the analyst doesn't have this

#I think we do this post-hoc and we can remove this.
# Extract the probe SNPs
${plink1.90} \
	--bfile ${bfile} \
	--extract resources/qc/probesnps.txt \
	--recode A \
	--out ${bfile}_probesnps


# Make GRMs and calculate PCs

gunzip -f resources/hapmap3_autosome.snplist.gz

${plink1.90} \
	--bfile ${bfile} \
	--extract resources/hapmap3_autosome.snplist \
	--maf 0.01 \
	--make-grm-bin \
	--out ${grmfile_all}

gzip -f resources/hapmap3_autosome.snplist

if [ "${family}" -eq "yes" ]
then
	Rscript resources/grm_relateds.R ${grmfile_all} ${grmfile_relateds} 0.125
elif [ "${family}" -eq "no" ]
then
	${plink1.90} \
		--grm-bin ${grmfile_all} \
		--rel-cutoff 0.125 \
		--make-grm-bin \
		--out ${grmfile_all}

	${gcta64} \
		--grm ${grmfile_all} \
		--pca 10 \
		--out ${genetic_pcs} \
		--thread-num ${nthreads}
else 
	echo "Error: Set family flag in config to yes or no"
	exit 1
fi



# Change SNVids to chr:position:{SNP/INDEL}

cp ${bfile}.bim ${bfile}.bim.original
awk '{if (length($5) == "1" && length($6) == "1") print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' <${bfile}.bim.original > ${bfile}.bim


# Generate meQTL format

#Generate plink.raw / plink.frq #freq files are required to determine effect allele
${plink1.90} --bfile ${bfile} --recode A --out ${bfile} --freq

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

./resources/make_meqtl_format.py ${bfile}



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

