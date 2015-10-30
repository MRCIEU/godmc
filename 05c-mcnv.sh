#!/bin/bash

# create batch scripts
# e.g to run:

# $ ./scriptname 1

set -e
source config

batch_number=${1}
re='^[0-9]+$'
if ! [[ $batch_number =~ $re ]] ; then
	echo "error: Batch variable is not a number"
	echo "Usage: ${0} [batch number]"
	exit 1
fi
exec &> >(tee ${mcnv_logfile}_${batch_number})

geno="${tabcnv}.tab.${batch_number}"
phen="${methylation_adjusted_pcs}.txt"
cov="NULL"
threshold="${soft_threshold}"
out="${matrixeqtl_mcnv_dir}/${resname}.${batch_number}.RData"


nbatch=(${tabcnv}.tab.*)
nbatch=${#nbatch[@]}

echo "Performing meQTL analysis on CNVs batch ${batch_number} of ${nbatch}"
Rscript resources/genetics/matrixeqtl.R ${geno} ${phen} ${cov} ${threshold} ${out}

echo "Successfully completed ${batch_number} of ${nbatch}"
