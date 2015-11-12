#!/bin/bash

set -e
source config
exec &> >(tee ${section_08_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else 
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${section_08_dir}/results \
		${section_08_dir}/qqplot
fi

echo "Successfully performed EWAS"
