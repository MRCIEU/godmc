#!/bin/bash

set -e
source config

# Check that 'successful' appears in the log files 

if grep -i -q "success" ${check_data_logfile}; then
	echo "01-check_data.sh completed successfully."
else
	echo "Problem: 01-check_data.sh did not complete successfully"
fi


if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -i -q "success" ${phenotype_data_logfile}; then
		echo "02a-phenotype_data.sh completed successfully."
	else
		echo "Problem: 02a-phenotype_data.sh did not complete successfully"
	fi
fi


if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -i -q "success" ${height_prediction_logfile}; then
		echo "02b-height_prediction.sh completed successfully."
	else
		echo "Problem: 02b-height_prediction.sh did not complete successfully"
	fi
fi


if grep -i -q "success" ${snp_data_logfile}; then
	echo "03a-snp_data.sh completed successfully."
else
	echo "Problem: 03a-snp_data.sh did not complete successfully"
fi


if grep -i -q "success" ${convert_snp_format_logfile}; then
	echo "03b-convert_snp_format.sh completed successfully."
else
	echo "Problem: 03b-convert_snp_format.sh did not complete successfully"
fi


if grep -i -q "success" ${methylation_variables_logfile}; then
	echo "04a-methylation_variables.sh completed successfully."
else
	echo "Problem: 04a-methylation_variables.sh did not complete successfully"
fi


if grep -i -q "success" ${methylation_adjustment1_logfile}*; then
	echo "04b-methylation_adjustment1.sh completed successfully."
else
	echo "Problem: 04b-methylation_adjustment1.sh did not complete successfully"
fi


if grep -i -q "success" ${methylation_pcs_logfile}; then
	echo "04c-methylation_pcs.sh completed successfully."
else
	echo "Problem: 04c-methylation_pcs.sh did not complete successfully"
fi


if grep -i -q "success" ${methylation_adjustment2_logfile}*; then
	echo "04d-methylation_adjustment2.sh completed successfully."
else
	echo "Problem: 04d-methylation_adjustment2.sh did not complete successfully"
fi

if grep -i -q "success" ${convert_methylation_format_logfile}*; then
	echo "04e-convert_methylation_format.sh completed successfully."
else
	echo "Problem: 04e-convert_methylation_format.sh did not complete successfully"
fi


nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -i "success" ${mqtl_logfile}* | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "05a-mqtl.sh completed successfully for all batches"
else
	echo "Problem: 05a-mqtl.sh only ${nsuccess} of ${nbatch} mQTL batches completed"
fi


nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -i "success" ${vmqtl_logfile}* | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "05b-vmqtl.sh completed successfully for all batches"
else
	echo "Problem: 05b-vmqtl.sh only ${nsuccess} of ${nbatch} variance-mQTL batches completed"
fi


nbatch=(${tabcnv}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -i "success" ${mcnv_logfile}* | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "05c-mcnv.sh completed successfully for all batches"
else
	echo "Problem: 05c-mcnv.sh only ${nsuccess} of ${nbatch} variance-mQTL batches completed"
fi


if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -i -q "success" ${ewas_logfile}; then
		echo "06-ewas.sh completed successfully."
	else
		echo "Problem: 06-ewas.sh did not complete successfully"
	fi
fi


if grep -i -q "success" ${gwas_aar_logfile}; then
	echo "07a-gwas_aar.sh completed successfully."
else
	echo "Problem: 07a-gwas_aar.sh did not complete successfully"
fi


if grep -i -q "success" ${gwas_smoking_logfile}; then
	echo "07b-gwas_smoking.sh completed successfully."
else
	echo "Problem: 07b-gwas_smoking.sh did not complete successfully"
fi


if grep -i -q "success" ${gwas_cellcount_entropy_logfile}; then
	echo "07c-gwas_cellcount_entropy.sh completed successfully."
else
	echo "Problem: 07c-gwas_cellcount_entropy.sh did not complete successfully"
fi


nsuccess=`grep -i "success" ${gwas_cellcounts_logfile}* | wc -l`
if [ "${nsuccess}" = "22" ]; then
	echo "07d-gwas_cellcounts.sh completed successfully on all chromosomes."
else
	echo "Problem: 07d-gwas_cellcounts.sh only complete for ${nsuccess} of 22 chromosomes"
fi

