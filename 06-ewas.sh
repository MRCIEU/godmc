#!/bin/bash

set -e
source config

if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else 
	R --no-save --args ${methylation_adjusted_pcs}.RData ${ewas_transformed} ${ewas_results} < resources/methylation/ewas.R
fi

