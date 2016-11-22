#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_04a_logfile})
print_version

# Predict smoking
echo "Estimating smoking"
Rscript resources/smoking/smoking_predictor.R \
	${betas} \
	${bfile}.fam \
	${smoking_pred} \
	${smoking_pred_plot} \
	${smoking_pred_SD} \
	${covariates} \
	${section_10_dir}/gwas_list.txt
	
# Estimate cell counts
echo "Processing cell counts for GWA"
if [ "${cellcounts_required}" == "NULL" ]
then
echo "Warning: Please change 'cellcounts_required' to yes or no in the config file."
exit 1
fi

if [ "${predicted_cellcounts}" == "NULL" ] && [ "${predicted_cellcounts_type}" != "NULL" ]
then
echo "Warning: Please change 'predicted_cellcounts_type' to NULL in the config file."
exit 1
fi

if [ "${cellcounts_required}" = "yes" ]
then
	if [ "${predicted_cellcounts}" == "NULL" ] && [ "${predicted_cellcounts_type}" == "NULL" ]
	then
		echo "Estimating cell counts using Houseman prediction"
		Rscript resources/cellcounts/estimate_cellcountsbybeta.R \
			${betas} \
			${cellcounts_gwa} \
			${cellcounts_reference} \
			${bfile}.fam
	fi
fi	

if [ -f "${predicted_cellcounts}" ] && [ "${predicted_cellcounts_type}" == "houseman" ]
then
	echo "Using the cellcounts provided in ${predicted_cellcounts} for GWA."
	cp ${predicted_cellcounts} ${cellcounts_gwa}
	
elif [ -f "${predicted_cellcounts}" ] && [ "${predicted_cellcounts_type}" != "houseman" ] && [ "${predicted_cellcounts_type}" != "NULL" ]
then
	echo "Your study will not be used for the cellcounts GWA meta-analysis - only Houseman predicted cellcounts are used"
else
	echo "Warning: The predicted cellcounts file ${predicted_cellcounts} doesn't exist. You have specified that cell counts are required. Please set 'predicted_cellcounts' to NULL if you want them to be estimated now, or specify a path to a file with the pre-specified cell counts."		
fi
	
if [ -f "${cellcounts_gwa}" ]
then
	echo "Transforming Houseman predicted cell counts for GWA"
	
        Rscript resources/genetics/create_cellcounts_plink.R \
		${cellcounts_gwa} \
		${cellcounts_plink_raw} \
		${intersect_ids_plink} \
		${cellcounts_plot} \
		${covariates} \
		${cellcounts_plink} \
		${cellcounts_SD} \
		${cellcounts_tf} \
		${cellcounts_entropy} \
		${smoking_pred} \
		${cellcounts_tf_smok} \
		${cellcounts_plink_smokadj}

elif [ "${cellcounts_required}" == "no" ]
then
	cellcounts_gwa="NULL"
else
	cellcounts_gwa="NULL"
	echo "'cellcounts_required' are set to yes but are not estimated with Houseman"
fi

# Estimate cell counts for use as covariates
echo "Processing cell counts for use as covariates"
if [ "${cellcounts_required}" = "yes" ]
then
	if [ "${predicted_cellcounts}" != "NULL" ] && [ "${measured_cellcounts}" == "NULL" ]
	then
		echo "Using the cellcounts provided in ${predicted_cellcounts} as covariates."
		cp ${predicted_cellcounts} ${cellcounts_cov}
	fi
	if [ -f "${measured_cellcounts}" ]
	then
		echo "Using the cellcounts provided in ${measured_cellcounts} as covariates."
		cp ${measured_cellcounts} ${cellcounts_cov}
	fi
    if [ "${predicted_cellcounts}" == "NULL" ] && [ "${predicted_cellcounts_type}" == "NULL" ]
	then
		echo "Using the cellcounts predicted from betas in ${cellcounts_gwa} as covariates."
		cp ${cellcounts_gwa} ${cellcounts_cov}
    fi
elif [ "${cellcounts_required}" = "no" ]
then
	cellcounts_cov="NULL"
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

# Organise covariates
echo "Organising covariates"
Rscript resources/genetics/covariates.R \
	${covariates} \
	${pcs_all} \
	${cellcounts_cov} \
	${smoking_pred}.txt \
	${bfile}.fam \
	${covariates_combined}

# GWAS Covariates
echo "Generating GWA covariates"
Rscript resources/genetics/create_covariates_files.R \
	${covariates_combined}.txt \
	${age_pred}.txt \
	${smoking_pred}.txt \
	${bfile}.fam \
	${gwas_covariates} \
	${covariates}

# GEMMA files

if [ -f "${cellcounts_gwa}" ]
then
echo "Preparing inputfiles for GEMMA"
Rscript resources/genetics/gemma_files.R \
	${grmfile_all} \
	${cellcounts_tf} \
	${section_12_dir}/cellcounts_columns.txt
fi

echo "Successfully created methylation-related variables"