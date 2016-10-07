#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16b_logfile})
print_version

# Get the control CPG
#positive_control_cpg="cg07959070"
#positive_control_snp_chr="22"
#positive_control_snp_pos="50053871"
#positive_control_snp_window="100000"
#positive_control_threshold="0.001"

colno=`awk -F' ' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' < ${methylation_processed_dir}/gcta.methylation.subset.1e-05_cg07[5-9].ge1.2.txt |grep -w ${positive_control_cpg} |awk '{print $1}'`

if [ ! "$colno" -gt "2" ];
then
	echo "The positive control CPG ${positive_control_cpg} appears to be missing. Please check."
	
fi

# Perform gcta

echo "Perform GWAS on ${positive_control_cpg}"

echo "CpG" "CHR" "SNP" "BP" "EA" "NEA" "EAF" "BETA" "SE" "P" | perl -pe 's/ /\t/g' > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control

for chr in `seq 1 23`;
${gcta} \
	--bfile ${bfile} \
	--mlma \
	--reml-no-constrain \
    --grm ${grmfile_all}_minus_chr${chrno} \
	--reml-no-constrain \
	--pheno ${methylation_processed_dir}/gcta.methylation.subset.1e-05_cg07[5-9].ge1.2.txt \
	--mpheno $colno
	--qcovar ${covariates_combined}.gcta.numeric \
	--covar ${covariates_combined}.gcta.factor \
	--out ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr \
	--chr $chr
	--thread-num ${nthreads}

tail -n +2 ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr >> ${section_16_dir}/gcta.${positive_control_cpg}.positive.control
rm ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr
fi


tr -s " " < ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.mlma | gzip -c > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.mlma.gz
rm ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.mlma

# make manhattan and qq plots

Rscript resources/genetics/plot_gwas.R \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold} \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control












