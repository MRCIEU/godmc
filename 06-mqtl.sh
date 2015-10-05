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
	exit 1
fi

geno="${tabfile}.tab.${batch_number}"

if [ "${unrelated}" -eq "no" ]
then
	phen="${methylation_rt_poly}"
else 
	phen="${methylation_rt}"
fi

cov="${covariates_combined}"
threshold=${soft_threshold}
out="${matrixeqtl_mqtl_dir}/${resname}.${batch_number}.RData"


R --no-save --args ${geno} ${phen} ${cov} ${threshold} ${out} < resources/genetics/matrixeqtl.R
