#!/bin/bash

set -e
source config


# Estimate cell counts
if [ "${cellcounts_required}" = "yes" ]
then
	R --no-save --args ${betas} ${cellcounts} < resources/cellcounts/estimate_cellcounts.R
elif [ "${cellcounts_required}" = "no" ]
then
	cellcounts="NULL"
else
	echo "'cellcounts_required' should be set to yes or no"
	exit 1
fi

# Adjust betas to rank transform and adjust for cell counts
R --no-save --args ${betas} ${cellcounts} ${methylation_rt} ${methylation_rt_cc} ${methylation_rt_cc_sq} ${nthreads} < resources/cellcounts/houseman.R


# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
if [ "${unrelated}" = "no" ]
then
	R --no-save --args ${betas} ${grmfile_relateds} ${methylation_rt_poly} ${nthreads} < resources/relateds/methylation_relateds.R
fi

