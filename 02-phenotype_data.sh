#!/bin/bash

set -e
source config
exec &> >(tee ${phenotype_data_logfile})

if [ ! "${phenotypes}" = "NULL" ]
then
Rscript	resources/phenotypes/phenotype_transform.R \
    ${intersect_ids_plink} \ 
    ${phenotypes} \ 
    ${ewasplot} \
    ${covariates} \
    ${ewastransformed} \
    ${ewas_sd} \
    ${ewas_plink} \
    ${ewas_file} \ 
fi

