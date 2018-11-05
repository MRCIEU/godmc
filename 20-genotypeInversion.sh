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

echo "Inferring inversions"

Rscript --vanilla resources/inversions/genotypeInversions.R ${inv_processed_dir}/SNPsonly ${inv_processed_dir} ${section_20_dir} ${nthreads}

echo "Successfully inferred inversions"
