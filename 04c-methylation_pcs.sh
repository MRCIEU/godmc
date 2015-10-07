#!/bin/bash

set -e
source config

# Calculate PCs for methylation, save in matrixeqtl format

R --no-save --args ${methylation_rt}.RData ${meth_pc_cutoff} ${covariates} ${meth_pcs} < resources/methylation/methylation_pcs.R


# Perform matrixeqtl on meth_pcs

R --no-save --args ${tabfile_hm3} ${meth_pcs} ${nongenetic_meth_pcs} ${meth_pc_threshold} ${nthreads} < resources/methylation/genetic_meth_pcs.R