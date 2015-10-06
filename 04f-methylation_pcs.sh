#!/bin/bash

set -e
source config

# Calculate PCs for methylation, save in matrixeqtl format

R --no-save --args ${methylation_rt}.RData ${n_meth_pcs} ${meth_pcs} < resources/methylation/methylation_pcs.R

# Perform matrixeqtl on meth_pcs

R --no-save --args ${tabfile_hm3} ${meth_pcs} ${nongenetic_meth_pcs} ${meth_pc_threshold} < resources/methylation/genetic_meth_pcs.R