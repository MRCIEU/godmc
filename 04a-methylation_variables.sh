#!/bin/bash

set -e
source config





# Estimate cell counts
if [ "${cellcounts_required}" = "yes" ]
then
	if [ "${provided_cellcounts}" = "NULL" ]
	then
		R --no-save --args ${betas} ${cellcounts} "${cellcount_reference}" ${bfile}.fam < resources/cellcounts/estimate_cellcounts_meffil.R

        elif [ -f "${provided_cellcounts}" ]
	then
		echo "Using the cellcounts provided in ${provided_cellcounts}."
		cp ${provided_cellcounts} ${cellcounts}
	else
		echo "Error: The file ${provided_cellcounts} doesn't exist. You have specified that cell counts are required. Please set 'provided_cellcounts' to NULL if you want them to be estimated now, or specify a path to a file with the pre-specified cell counts."
	fi
	#R --no-save --args ${cellcounts} ${cellcounts_plink} ${home_directory}/processed_data/cellcounts/  < resources/genetics/create_cellcounts_plink.R
	R --no-save --args ${cellcounts} ${cellcounts_plink_raw} ${intersect_ids_plink} ${cellcounts_plot} ${covariates} ${cellcounts_plink} ${cellcounts_SD} < resources/genetics/create_cellcounts_plink.R ${home_directory}/processed_data/cellcounts/cellcounts_transform.Rout 2>&1

elif [ "${cellcounts_required}" = "no" ]
then
	cellcounts="NULL"
else
	echo "'cellcounts_required' should be set to yes or no"
	exit 1
fi

# Estimate age accelerated residuals
R --no-save --args ${betas} ${covariates} ${bfile}.fam ${age_pred} ${age_pred_plot} ${age_pred_SD} < resources/dnamage/dnamage.R

# Predict smoking
R --no-save --args ${betas} ${bfile}.fam ${smoking_pred} ${smoking_pred_plot} ${smoking_pred_SD} ${covariates} < resources/smoking/smoking_predictor.R

# Organise covariates
R --no-save --args ${covariates} ${pcs_all} ${cellcounts} ${smoking_pred}.txt ${bfile}.fam ${covariates_combined} < resources/genetics/covariates.R

# GWAS Covariates
R --no-save --args ${covariates_combined}.txt ${age_pred}.txt ${smoking_pred}.txt ${bfile}.fam
 ${gwas_covariates} < resources/genetics/create_covariate_files.R
