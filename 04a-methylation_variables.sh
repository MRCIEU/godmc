#!/bin/bash

set -e
source config
exec &> >(tee ${methylation_variables_logfile})


# Estimate cell counts
if [ "${cellcounts_required}" = "yes" ]
then
	if [ "${provided_cellcounts}" = "NULL" ]
	then
		echo "Estimating cell counts"
		Rscript resources/cellcounts/estimate_cellcountsbybeta.R \
			${betas} \
			${cellcounts} \
			${cellcounts_reference} \
			${bfile}.fam
		elif [ -f "${provided_cellcounts}" ]
	then
		echo "Using the cellcounts provided in ${provided_cellcounts}."
		cp ${provided_cellcounts} ${cellcounts}
	else
		echo "Error: The file ${provided_cellcounts} doesn't exist. You have specified that cell counts are required. Please set 'provided_cellcounts' to NULL if you want them to be estimated now, or specify a path to a file with the pre-specified cell counts."
	fi
	echo "Transforming cell counts"
	Rscript resources/genetics/create_cellcounts_plink.R \
		${cellcounts} \
		${cellcounts_plink_raw} \
		${intersect_ids_plink} \
		${cellcounts_plot} \
		${covariates} \
		${cellcounts_plink} \
		${cellcounts_SD} \
		${cellcounts_tf} \
		${cellcounts_entropy}

elif [ "${cellcounts_required}" = "no" ]
then
	cellcounts="NULL"
else
	echo "'cellcounts_required' should be set to yes or no"
	exit 1
fi

# Estimate age accelerated residuals
echo "Estimating age"
Rscript resources/dnamage/dnamage.R \
	${betas} \
	${covariates} \
	${bfile}.fam \
	${age_pred} \
	${age_pred_plot} \
	${age_pred_SD}

# Predict smoking
echo "Estimating smoking"
Rscript resources/smoking/smoking_predictor.R \
	${betas} \
	${bfile}.fam \
	${smoking_pred} \
	${smoking_pred_plot} \
	${smoking_pred_SD} \
	${covariates}

# Organise covariates
echo "Organising covariates"
Rscript resources/genetics/covariates.R \
	${covariates} \
	${pcs_all} \
	${cellcounts} \
	${smoking_pred}.txt \
	${bfile}.fam \
	${covariates_combined}

# GWAS Covariates
Rscript resources/genetics/create_covariates_files.R \
	${covariates_combined}.txt \
	${age_pred}.txt \
	${smoking_pred}.txt \
	${bfile}.fam \
	${gwas_covariates}

# GEMMA files
Rscript resources/genetics/gemma_files.R \
	${grmfile_all} \
	${cellcounts_tf} \
	${gwas_cellcounts_dir}/cellcounts_columns.txt \
	${gwas_covariates}.cellcounts

echo "Successfully created methylation-related variables"
