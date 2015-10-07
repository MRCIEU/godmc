#!/bin/bash

set -e
source config

batch_number=${1}



if [ -n "${1}" ]
then
	re='^[0-9]+$'
	if ! [[ $batch_number =~ $re ]] ; then
		echo "error: Batch variable is not a number"
		exit 1
	fi
	i=${1}
	echo "Running batch ${i} of ${meth_chunks}"
else
	i="NA"
	echo "Running entire set on a single node using ${nthreads} threads."
fi


if [ "${unrelated}" = "yes" ]
then
	echo "You have specified that the data is not family data. Adjusting only for covariates..."
	R --no-save --args ${betas} ${covariates_combined}.txt ${methylation_adjusted} ${nthreads} ${meth_chunks} ${i} < resources/methylation/adjust_covs.R
elif [ "${unrelated}" = "no" ]
then
	# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
	echo "You have specified that the data is family data. Adjusting for pedigree and covariates..."
	R --no-save --args ${betas} ${grmfile_relateds} ${covariates_combined}.txt ${methylation_adjusted} ${nthreads} ${meth_chunks} ${i} < resources/methylation/adjust_pedigree.R
fi
