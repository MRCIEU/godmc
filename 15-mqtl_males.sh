#!/bin/bash

set -e
source ./config

batch_number=${1}
re='^[0-9]+$'
if ! [[ $batch_number =~ $re ]] ; then
	echo "error: Batch variable is not a number"
	echo "Usage: ${0} [batch number]"
	exit 1
fi
exec &> >(tee ${section_15_logfile}${batch_number})
print_version

if [ ! -f "${methylation_adjusted_pcs}_males.txt" ]
then
	echo "No male-only data. Exiting Successfully"
	exit
fi

geno="${tabfile}.tab.${batch_number}"
phen="${methylation_adjusted_pcs}_males.txt"
cov="NULL"
threshold="0.01"
out="${section_15_dir}/${resname}.${batch_number}.RData"


nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}

echo "Performing meQTL analysis batch ${batch_number} of ${nbatch}"
Rscript resources/genetics/matrixeqtl.R ${geno} ${phen} ${cov} ${threshold} ${out}

echo "Successfully completed ${batch_number} of ${nbatch}"
