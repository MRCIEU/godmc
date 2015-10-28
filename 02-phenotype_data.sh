#!/bin/bash

set -e
source config


if [ ! "${phenotypes}" = "NULL" ]
then
	R --no-save --args ${intersect_ids_plink} ${phenotypes} ${ewasplot} ${covariates} ${ewastransformed} ${ewas_sd} ${ewas_plink} ${ewas_file} < resources/phenotypes/phenotype_transform.R 
fi
