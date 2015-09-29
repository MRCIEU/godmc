#!/bin/bash

set -e
source config




# Use output from meffil to normalise data and generate betas in RData format

# Adjust betas to rank transform and adjust for cell counts
R --no-save --args ${betas} ${methylation_rt} ${methylation_rt_cc} ${methylation_rt_cc_sq} ${cellcounts} ${nthreads} < resources/cellcounts/houseman.R

# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
if [ "${family}" -eq "yes" ]
then
	R --no-save --args ${betas} ${grmfile_relateds} ${methylation_rt_poly} ${nthreads} < resources/relateds/methylation_relateds.R
fi

# Organise covariates
R --no-save --args ${covariates} ${pcs_all} ${cellcounts} ${covariates_combined} < resources/genetics/covariates.R
