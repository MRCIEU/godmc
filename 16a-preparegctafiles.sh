#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16a_logfile})
print_version


# Extract the SNP-CpG lists

	sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
mget cis_trans.1e-05_*.ge${no}.allcohorts.txt.gz
mget cis_trans.1e-05_*.ge${no}.allcohorts.probes
EOF

	mv cis_trans.1e-05_*.ge${no}.allcohorts.txt.gz ${home_directory}/resources/phase2
	mv cis_trans.1e-05_*.ge${no}.allcohorts.probes ${home_directory}/resources/phase2


# Create beta matrix for Plink

echo "Creating files for plink"

# Subset the adjusted beta matrix from 04d per probe set 
Rscript resources/phase2/extract_relevant_probes.R \
	${methylation_adjusted_pcs}.RData \
	${bfile}.fam \
	${phase2_assoclist} \
	${phase2_betas}

cut -d " " -f 1-2 ${phase2_betas}1 | sed 1d > keeplist.txt

norig=`cat ${bfile}.fam | wc -l`
nnow=`cat keeplist.txt | wc -l`

echo "${norig} samples in original genotype file"
echo "${nnow} samples with both genotype and methylation data"

${plink} \
	--bfile ${bfile} \
	--keep keeplist.txt \
	--maf 0.01 \
	--make-bed \
	--out ${bfile}_phase2


echo "Calculating MAF for individuals with CpG data"

${plink} \
    --bfile ${bfile}_phase2 \
    --freq gz \
    --out ${section_16_dir}/data

zcat ${section_16_dir}/data.frq.gz | awk '{ print $1,$2,$3,$4,$5,$6/2 }' | sed 1d | gzip -c > temp.gz
zcat temp.gz | awk '{ print $2 }' | gzip > ${section_16_dir}/snplist.txt.gz
mv temp.gz ${section_16_dir}/data.frq.gz
rm keeplist.txt


echo "Successfully completed script 16a"

