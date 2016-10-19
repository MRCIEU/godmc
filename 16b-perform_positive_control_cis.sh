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

colno=$((colno - 2))

#find index cisSNP
echo "Find cisindexSNP for ${positive_control_cpg}"

${gcta} \
	--bfile ${bfile} \
	--mlma \
	--grm ${grmfile_all}_minus_chr$positive_control_snp_chr \
	--pheno ${methylation_processed_dir}/gcta.methylation.subset.1e-05_cg07[5-9].ge1.2.txt \
	--mpheno $colno \
	--qcovar ${covariates_combined}.gcta.numeric \
	--covar ${covariates_combined}.gcta.factor \
	--chr $positive_control_snp_chr \
	--out ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr \
	--thread-num ${nthreads}

minp=`awk '{print $9}' <${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.mlma |sort -g |tail -2 |head -1`
cissnp=`awk '$9=='$minp'{print $2}' <${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cissnp`

echo "extract cissnp"

${plink} \
	--bfile ${bfile} \
	--allow-no-sex \
	--extract ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cissnp \
	--recode A \
	--out ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr

echo "include cisSNP in covariate file"

awk 'NR>1 {print $7}' < ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.raw > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.tmp
paste ${covariates_combined}.gcta.numeric ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.tmp |perl -pe 's/\t/ /g' > ${covariates_combined}.gcta.numeric.cis

# Perform gcta

echo "Perform GWAS on ${positive_control_cpg}"

echo "Chr" "SNP" "bp" "A1" "A2" "Freq" "b" "se"	"p" |perl -pe 's/ /\t/g' > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis

for chr in `seq 1 22`;
do

echo "chromosome $chr"

${gcta} \
	--bfile ${bfile} \
	--mlma \
	--grm ${grmfile_all}_minus_chr$chr \
	--pheno ${methylation_processed_dir}/gcta.methylation.subset.1e-05_cg07[5-9].ge1.2.txt \
	--mpheno $colno \
	--qcovar ${covariates_combined}.gcta.numeric.cis \
	--covar ${covariates_combined}.gcta.factor \
	--chr $chr \
	--out ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr.cis \
	--thread-num ${nthreads}

tail -n +2 ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr.cis.mlma >> ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis

done

tr -s " " < ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis | gzip -c > ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis.gz

# make manhattan and qq plots

Rscript resources/genetics/plot_gwas.R \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis.gz \
	9 \
	1 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold} \
	${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis



if [ -f "${section_16_dir}/gcta.${positive_control_cpg}.positive.control.gz" ]
then 
rm ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr*.mlma
fi

echo "Successfully completed script 16b"






