#!/bin/bash

set -e
source ./config

exec &> >(tee ${section_21_logfile})
print_version

geno="${inv_processed_dir}/inversionQTL.txt"
phen="${methylation_adjusted_pcs}.txt"
out="${section_21_dir}/invmeqtl.txt.gz"

echo "Performing inversion meQTL analysis"
Rscript resources/inversions/run_analysis_inversions.R ${geno} ${phen} ${out}

echo "Successfully completed inversion meQTL analysis"
