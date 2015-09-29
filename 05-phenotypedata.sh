#!/bin/bash

set -e
source config


# Use output from meffil to normalise data and generate betas in RData format
if [ -z "${provided_cellcounts}" ]
then
	R --no-save --args 
else 
	cp ${provided_cellcounts} ${cellcounts}
fi

# Predict smoking
R --no-save --args ${betas} ${smoking_pred} < resources/smoking/smoking_predictor.R

# Organise covariates
# ADD SMOKING PREDICTION TO COVARIATES
R --no-save --args ${covariates} ${pcs_all} ${cellcounts} ${smoking_pred} ${covariates_combined} < resources/genetics/covariates.R

# Estimate age accelerated residuals
R --no-save --args ${betas} ${covariates} ${age_pred} < resources/dnamage/dnamage.R

# Height and BMI adjustments
