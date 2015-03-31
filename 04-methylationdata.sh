#!/bin/bash

set -e
source config


# Use output from meffil to normalise data and generate betas in RData format


# Adjust betas to rank transform and adjust for cell counts
Rscript resources/houseman/houseman.R ${beta} ${rnbeta} ${ccrnbeta} ${cellcounts} ${nthreads}

# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
if [ "${family}" -eq "yes" ]
then
	Rscript resources/methylation_relateds.R ${ccrnbeta} ${grmfile_all} ${ccrnfambeta} 0.05
fi

# Use R-GADA to generate structural variant data from methylation intensities


# Filter methylation data on Naeem et al list (categories) and non-variant probes (threshold) and non-unique positions


save
 structural variants in processed_data
 normalised methylation in processed_data
 
