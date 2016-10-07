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

${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--reml-no-constrain \
	--grm ${grmfile_all} \
	--pheno ${methylation_processed_dir}/gcta.methylation.subset.1e-05_cg07[5-9].ge1.2.txt \
	--mpheno $colno
	--qcovar ${covariates_combined}.gcta.numeric \
	--covar ${covariates_combined}.gcta.factor \
	--out ${section_16_dir}/gcta.${positive_control_cpg}.positive.control \
	--thread-num ${nthreads}




tr -s " " < ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.loco.mlma | gzip -c > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.loco.mlma.gz
rm ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.loco.mlma

# make manhattan and qq plots

Rscript resources/genetics/plot_gwas.R \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control.loco.mlma.gz \
	9 \
	1 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold} \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control











	
