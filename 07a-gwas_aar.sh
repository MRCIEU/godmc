#!/bin/bash

set -e
source config

# Setup covariates
# One covariate - sex

${gcta} --bfile ${bfile} --mlma-loco --pheno ${age_pred}.aar.plink --qcovar ${gwas_covariates}.aar --out ${gwas_aar_dir}/gwas

# Setup phenotype
# Run GCTA

