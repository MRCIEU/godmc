#!/bin/bash

set -e

# Provision for having dosage data and convert to best guess if the analyst doesn't have this


source config


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
	Rscript resources/grm_relateds.R ${grmfile_all} ${grmfile_relateds} 0.05
elif [ "${family}" -eq "no" ]
then
	${plink1.90} \
		--grm-bin ${grmfile_all} \
		--rel-cutoff 0.05 \
		--make-grm-bin \
		--out ${grmfile_all}

	${gcta64} \
		--grm ${grmfile_all} \
		--pca 10 \
		--out ${genetic_pcs} \
		--thread-num ${nthreads}
fi



# Change SNPs to chr:position

cp ${bfile}.bim ${bfile}.bim.original
awk ${bfile}.bim.original '{ print $1, "chr"$1":"${4}, $3, $4, $5, $6 }' > ${bfile}.bim


# Generate meQTL format

${plink1.90} --bfile ${bfile} --recode A --out ${bfile}
./resources/make_meqtl_format.py ${bfile}



# make GRM from hm3 SNPs
# if family data then make pedigree matrix
# if unrelateds then calculate principal components
# convert bim file to change rs IDs to chr<x>:<pos>
# remove cryptically related individuals

# create matrix eqtl format 
## hash has this already use --recode option


