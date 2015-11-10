#!/bin/bash

set -e
source config
exec &> >(tee ${section_08_logfile})

if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else 
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${ewas_results_dir}/results \
		${ewas_results_dir}/qqplot
fi

echo "Successfully performed EWAS"
