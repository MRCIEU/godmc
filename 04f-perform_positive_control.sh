#!/bin/bash

set -e
source config
exec &> >(tee ${section_04f_logfile})
print_version

# Get the control CPG
awk -v cpg=${positive_control_cpg} '{ if(NR == 1 || $1 == cpg) { print $0 }}' ${methylation_adjusted_pcs}.txt > ${methylation_adjusted_pcs}.positive_control

nrow=`cat ${methylation_adjusted_pcs}.positive_control | wc -l`
if [ "$nrow" -lt "2" ];
then
	echo "The positive control CPG ${positive_control_cpg} appears to be missing. Please check."
	exit
fi

Rscript resources/genetics/make_control.R \
	${methylation_adjusted_pcs}.positive_control \
	${intersect_ids_plink} \
	${methylation_adjusted_pcs}.positive_control.plink

# Perform plink

${plink} \
	--bfile ${bfile} \
	--pheno ${methylation_adjusted_pcs}.positive_control.plink \
	--assoc \
	--out ${section_04_dir}/positive_control_${positive_control_cpg}

tr -s " " < ${section_04_dir}/positive_control_${positive_control_cpg}.qassoc | gzip -c > ${section_04_dir}/positive_control_${positive_control_cpg}.qassoc.gz
rm ${section_04_dir}/positive_control_${positive_control_cpg}.qassoc

# make manhattan and qq plots

Rscript resources/genetics/plot_gwas.R \
	${section_04_dir}/positive_control_${positive_control_cpg}.qassoc.gz \
	9 \
	1 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold} \
	${section_04_dir}/positive_control_${positive_control_cpg}

echo "Successfully completed check"
