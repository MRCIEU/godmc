#!/bin/bash

set -e
source config
if [ "${EWAS_phenotypes}" = "yes" ]
then
R --no-save --args ${intersect_ids} ${EWASphenotypes} ${EWASplot} ${covariates} ${EWAStransformed} ${EWAS_SD} < resources/phenotypes/phenotype_transform.R > ${home_directory}/processed_data/phenotypes/phenotype_transform.Rout 2>&1
fi
# Height and BMI adjustments
