#!/bin/bash

set -e
source config
if [ "${EWAS_phenotypes}" = "yes" ]
then
R --no-save --args ${EWASphenotypes} ${home_directory}/processed_data/phenotypes/EWAS.phenotypes.pdf ${covariates} ${home_directory}/processed_data/phenotypes/EWAS.phenotypes.transformed.txt 5 < resources/phenotypes/phenotype_transform.R > ${home_directory}/processed_data/phenotypes/phenotype_transform.Rout 2>&1
fi
# Height and BMI adjustments
