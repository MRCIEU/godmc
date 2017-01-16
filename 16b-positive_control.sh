#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16b_logfile})
print_version

# Get the control CPG - in parameters file
#positive_control_cpg="cg07959070"
#positive_control_snp_chr="22"
#positive_control_snp_pos="50053871"
#positive_control_snp_window="100000"
#positive_control_threshold="0.001"

# subset of positive control
pos_sub ="cg07[5-9]"

# 128 subsets
colno=`awk -F' ' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' < ${methylation_processed_dir}/plink.methylation.subset.1e-05_${pos_sub}.ge1.2.txt |grep -w ${positive_control_cpg} |awk '{print $1}'`

if [ ! "$colno" -gt "2" ];
then
	echo "The positive control CPG ${positive_control_cpg} appears to be missing. Please check."
	
fi

colno=$((colno - 2))

#find index cisSNP
echo "Find cisindexSNP for ${positive_control_cpg}"

#i="1e-05"
#no="1.2"

# GWAS of positive control probe
${plink} \
    --bfile ${bfile} \
    --extract ${genetic_processed_dir}/cis_trans.1e-05_${pos_sub}.ge1.2.allcohorts.snps \
    --pheno ${methylation_processed_dir}/plink.methylation.subset.1e-05_${pos_snp}.ge1.2.txt \
    --mpheno $colno \
    --out ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr \
    --allow-no-sex \
    --assoc \
    --threads ${nthreads}

# Plink ouput: Chr   SNP  bp     A1         A2       Freq b       se     p

# Extract most strongly associated cis SNP to positive control
minp=`awk '{print $9}' <${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.qassoc |sort -g |tail -2 |head -1`
cissnp=`awk '$9=='$minp'{print $2}' <${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.qassoc > ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cissnp`

echo "extract cissnp"

${plink} \
	--bfile ${bfile} \
	--allow-no-sex \
	--extract ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cissnp \
	--recode A \
	--out ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr

# Include cisSNP as a covariate in reg
	
	echo "include cisSNP in GWAS of positive control"

# Gen covar file? --covar
# Condition on SNP? --condition
	
######### COME BACK TO THIS	- written for covar files of GCTA
awk 'NR>1 {print $7}' < ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.raw > ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cissnp.tmp
paste ${covariates_combined}.plink ${section_16_dir}/plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.tmp |perl -pe 's/\t/ /g' > ${covariates_combined}.plink.cis

echo "Chr" "SNP" "bp" "A1" "A2" "Freq" "b" "se"	"p" |perl -pe 's/ /\t/g' > ${section_16_dir}/plink.${positive_control_cpg}.positive.control.cis
########

# Perform GWAS on positive control

echo "Perform GWAS on ${positive_control_cpg}"

${plink} \
    --bfile ${bfile} \
    --extract ${genetic_processed_dir}/cis_trans.1e-05_${pos.sub}.ge1.2.allcohorts.snps \
    --pheno ${methylation_processed_dir}/plink.methylation.subset.1e-05_${pos_snp}.ge1.2.txt \
    --mpheno $colno \
	--
    --out plink.${positive_control_cpg}.positive.control.chr$positive_control_snp_chr.cis \
    --allow-no-sex \
    --assoc \
    --threads ${nthreads}


# Format Results

############# Amend this - from 16c-run_mqtl.sh
sed 's/^/'$probe'/' <${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
            tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${section_16_dir}/plink.${i}\_${j}.ge${no}.txt
        
            rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
            rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc
            rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.log
            rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps	

tail -n +2 ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr$chr.cis.mlma >> ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.cis
###############

# Compress
tr -s " " < ${section_16_dir}/plink.${positive_control_cpg}.positive.control.cis | gzip -c > ${section_16_dir}/plink.${positive_control_cpg}.positive.control.cis.gz

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
	${section_16_dir}/plink.${positive_control_cpg}.positive.control.cis


######### I think this is here to remove non cis effect positive control? (GCTA)
if [ -f "${section_16_dir}/plink.${positive_control_cpg}.positive.control.gz" ]
then 
rm ${section_16_dir}/gcta.${positive_control_cpg}.positive.control.chr*.mlma
fi
#########

echo "Successfully completed script 16b"






