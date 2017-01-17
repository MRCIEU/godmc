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
Rscript ./resources/phase2/extract_relevant_probes.R \
	${methylation_adjusted_pcs}.RData \
	${bfile}.fam \
	${phase2_assoclist} \
	${phase2_betas}

echo "Calculating MAF"
${plink} \
    --bfile ${bfile} \
    --freq gz \
    --out ${section_16_dir}/data

zcat ${section_16_dir}/data.frq.gz | awk '{ print $1,$2,$3,$4,$5,$6 }' | gzip -c > temp
mv temp ${section_16_dir}/data.frq.gz


echo "Successfully completed script 16a"

