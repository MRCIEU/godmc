#!/bin/bash

set -e
source config
exec &> >(tee ${check_results_logfile})

# Check that 'successful' appears in the log files 

echo ""
echo "Checking log files"
echo "If any log files are flagged as having problems please re-run the relevant scripts"
echo ""

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
	echo "Problem: 05c-mcnv.sh only ${nsuccess} of ${nbatch} CNV-mQTL batches completed"
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



## Check results files

echo ""
echo "Checking results files"
echo "If any results files are flagged as having problems please re-run the relevant scripts"
echo ""


# EWAS
if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -q "Height" ${ewas_results}.phenotype_list; then
		if grep -l -q "Height" ${ewas_results}*; then
			echo "Height EWAS results present"
		else
			echo "Problem: Height EWAS results absent"
		fi
	fi
	if grep -q "BMI" ${ewas_results}.phenotype_list; then
		if grep -l -q "BMI" ${ewas_results}*; then
			echo "BMI EWAS results present"
		else
			echo "Problem: BMI EWAS results absent"
			exit 1
		fi
	fi
fi


# smoking GWAS
if grep -l -q "loco" ${gwas_smoking_dir}/*; then
	echo "Smoking GWAS results present"
else
	echo "Problem: Smoking GWAS results absent"
fi

# aar GWAS
if grep -l -q "loco" ${gwas_aar_dir}/*; then
	echo "AAR GWAS results present"
else
	echo "Problem: AAR GWAS results absent"
fi

# cellcount_entropy GWAS
if grep -l -q "loco" ${gwas_cellcounts_dir}/*; then
	echo "Cellcount entropy GWAS results present"
else
	echo "Problem: Cellcount entropy GWAS results absent"
fi

# cellcounts GWAS
if [ -f "${gwas_cellcounts_dir}/cellcounts_columns.txt" ]; then
	echo "Cellcount annotation present"
else
	echo "Problem: Cellcount annotation is absent"
fi

nsuccess=`grep -l "assoc" ${gwas_cellcounts_dir}/* | wc -l`
if [ "${nsuccess}" = "22" ]; then
	echo "Cellcounts results present for all chromosomes"
else
	echo "Problem: Only ${nsuccess} of 22 cellcounts result files are present"
fi

# mqtl
nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -l "res" ${matrixeqtl_mqtl_dir}/* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} mQTL result files are present"
fi

if [ -f "${allele_ref}" ]; then
	echo "Allele reference file present"
else
	echo "Problem: Allele reference file is absent"
fi


# vmqtl
nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -l "res" ${matrixeqtl_vmqtl_dir}/* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All variance-mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} variance-mQTL result files are present"
fi

# mcnv
nbatch=(${tabcnv}.tab.*)
nbatch=${#nbatch[@]}
nsuccess=`grep -l "res" ${matrixeqtl_mcnv_dir}/* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All CNV-mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} CNV-mQTL result files are present"
fi

