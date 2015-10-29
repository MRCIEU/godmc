#!/bin/bash

set -e
source config
exec &> >(tee ${methylation_pcs_logfile})

# Calculate PCs for methylation, save in matrixeqtl format

echo "Calculating methylation PCs"
Rscript resources/methylation/methylation_pcs.R \
	${methylation_rt}.RData \
	${meth_pc_cutoff} \
	${normalised_phenotypes} \
	${meth_pcs}


# Perform matrixeqtl on meth_pcs

Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile_hm3} \
	${meth_pcs} \
	${nongenetic_meth_pcs} \
	${meth_pc_threshold} \
	${nthreads}