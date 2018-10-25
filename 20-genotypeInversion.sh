#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_20_logfile})
print_version

# Change annotation ####
#'#################################################################################

echo "Adapting plink files"


## Remove indels
cut -f2 ${bfile}.bim | grep "INDEL" > ${inv_processed_dir}/indels.txt
${plink} -bfile  ${bfile} --exclude ${inv_processed_dir}/indels.txt --make-bed --out ${inv_processed_dir}/SNPsonly


## Change SNP names to chr:pos
awk '{print $2, $1":"$4}' ${inv_processed_dir}/SNPsonly.bim > ${inv_processed_dir}/newAnnot.tab
${plink} -bfile  ${bfile} --update-name ${inv_processed_dir}/newAnnot.tab --make-bed --out ${inv_processed_dir}/newAnnot

## Select SNPs in list
${plink} -bfile ${inv_processed_dir}/newAnnot --extract ${inv_snps} --snps-only --make-bed --out ${inv_processed_dir}/filtSNPs

### Update annotation
${plink} -bfile ${inv_processed_dir}/filtSNPs --update-name $inv_snps --make-bed --out ${inv_processed_dir}/rsAnnot

echo "Inferring inversions"

Rscript --vanilla resources/inversions/genotypeInversions.R ${inv_processed_dir}/rsAnnot ${inv_processed_dir}

echo "Successfully inferred inversions"
