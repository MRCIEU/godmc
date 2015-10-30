#!/bin/bash

set -e
source config
exec &> >(tee ${methylation_pcs_logfile})

# Calculate PCs for methylation, save in matrixeqtl format

if [ ! "${phenotypes}" = "NULL" ]
then
	phenfile="${ewastransformed}"
else
	phenfile="NULL"
fi

echo "Calculating methylation PCs"
Rscript resources/methylation/methylation_pcs.R \
	${methylation_adjusted}.RData \
	${meth_pc_cutoff} \
	${phenfile} \
	${meth_pcs}


echo "Performing genetic analysis of methylation PCs"
Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile_hm3} \
	${meth_pcs} \
	${nongenetic_meth_pcs} \
	${meth_pc_threshold} \
	${nthreads}
