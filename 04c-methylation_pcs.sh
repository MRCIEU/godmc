#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_04c_logfile})
print_version

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
	${n_meth_pcs} \
	${phenfile} \
	${meth_pcs}


echo "Performing genetic analysis of methylation PCs"
Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile} \
	${meth_pcs} \
	${nongenetic_meth_pcs} \
	${meth_pc_threshold} \
	${nthreads}

echo "Successfully generated non-genetic methylation PCs"