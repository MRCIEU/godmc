#!/bin/bash

set -e
source config
exec &> >(tee ${ewas_logfile})

if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else 
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${ewas_results}
fi

echo "Successfully performed EWAS"
