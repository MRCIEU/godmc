#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16a_logfile})
print_version


# Extract the SNP-CpG lists

echo "Downloading list of putative associations"

sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
get phase2_list.tar
get phase2_list.tar.md5sum
EOF

echo "Checking download integrity"

md5sum -c phase2_list.tar.md5sum

echo "Extracting"

tar xf phase2_list.tar -C ${phase2_assoclist}
rm phase2_list.tar*


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

