#!/bin/bash

set -e
source config


checkFirstArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid section identifier"
	echo $"Need to specify a value from 01 to $e"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}


checkSecondArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $2 is not a valid action"
	echo $"Specify either 'check' or 'upload'"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}


check_logs_01 () {

	if grep -i -q "success" ${check_data_logfile}; then
		echo "01-check_data.sh completed successfully."
	else
		echo "Problem: 01-check_data.sh did not complete successfully"
		exit 1
	fi

}



check_logs_02 () {


	if grep -i -q "success" ${snp_data_logfile}; then
		echo "02a-snp_data.sh completed successfully."
	else
		echo "Problem: 02a-snp_data.sh did not complete successfully"
		echo "Please ensure previous steps have been performed and rerun 02a-snp_data.sh"
		exit 1
	fi


	if grep -i -q "success" ${convert_snp_format_logfile}; then
		echo "02b-check_snp_format.sh completed successfully."
	else
		echo "Problem: 02b-check_snp_format.sh did not complete successfully"
		exit 1
	fi

}



check_logs_03 () {

	if [ "${phenotypes}" = "NULL" ]
	then
		echo "No phenotypes specified"
		exit 1
	fi

	if grep -i -q "success" ${phenotype_data_logfile}; then
		echo "03a-phenotype_data.sh completed successfully."
	else
		echo "Problem: 03a-phenotype_data.sh did not complete successfully"
		exit 1
	fi

	if grep -i -q "success" ${height_prediction_logfile}; then
		echo "03b-height_prediction.sh completed successfully."
	else
		echo "Problem: 03b-height_prediction.sh did not complete successfully"
		exit 1
	fi

}


check_logs_04 () {

	if grep -i -q "success" ${methylation_variables_logfile}; then
		echo "04a-methylation_variables.sh completed successfully."
	else
		echo "Problem: 04a-methylation_variables.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${methylation_adjustment1_logfile}*; then
		echo "04b-methylation_adjustment1.sh completed successfully."
	else
		echo "Problem: 04b-methylation_adjustment1.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${methylation_pcs_logfile}; then
		echo "04c-methylation_pcs.sh completed successfully."
	else
		echo "Problem: 04c-methylation_pcs.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${methylation_adjustment2_logfile}*; then
		echo "04d-methylation_adjustment2.sh completed successfully."
	else
		echo "Problem: 04d-methylation_adjustment2.sh did not complete successfully"
		exit 1
	fi

	if grep -i -q "success" ${convert_methylation_format_logfile}*; then
		echo "04e-convert_methylation_format.sh completed successfully."
	else
		echo "Problem: 04e-convert_methylation_format.sh did not complete successfully"
		exit 1
	fi

}


check_logs_05 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`tail ${mqtl_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "05a-mqtl.sh completed successfully for all batches"
	else
		echo "Problem: 05a-mqtl.sh only ${nsuccess} of ${nbatch} mQTL batches completed"
		exit 1
	fi


	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`tail ${vmqtl_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "05b-vmqtl.sh completed successfully for all batches"
	else
		echo "Problem: 05b-vmqtl.sh only ${nsuccess} of ${nbatch} variance-mQTL batches completed"
		exit 1
	fi


	nbatch=`ls -l ${tabcnv}.tab.* | wc -l`
	nsuccess=`tail ${mcnv_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "05c-mcnv.sh completed successfully for all batches"
	else
		echo "Problem: 05c-mcnv.sh only ${nsuccess} of ${nbatch} CNV-mQTL batches completed"
		exit 1
	fi

}


check_logs_06 () {

	if [ "${phenotypes}" = "NULL" ]
	then
		echo "No phenotypes specified"
		exit 1
	fi

	if grep -i -q "success" ${ewas_logfile}; then
		echo "06-ewas.sh completed successfully."
	else
		echo "Problem: 06-ewas.sh did not complete successfully"
		exit 1
	fi

}



check_logs_07 () {

	if grep -i -q "success" ${gwas_aar_logfile}; then
		echo "07a-gwas_aar.sh completed successfully."
	else
		echo "Problem: 07a-gwas_aar.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${gwas_smoking_logfile}; then
		echo "07b-gwas_smoking.sh completed successfully."
	else
		echo "Problem: 07b-gwas_smoking.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${gwas_cellcount_entropy_logfile}; then
		echo "07c-gwas_cellcount_entropy.sh completed successfully."
	else
		echo "Problem: 07c-gwas_cellcount_entropy.sh did not complete successfully"
		exit 1
	fi

	nbatch=`wc -l ${gwas_cellcounts_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
	nsuccess=`tail ${gwas_cellcounts_logfile}_* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "07d-gwas_cellcounts.sh completed successfully for all cell types"
	else
		echo "problem: 07d-gwas_cellcounts.sh only complete for ${nsuccess} of ${nbatch} cell types"
		exit 1
	fi

	nsuccess=`grep -i "success" ${gwas_cellcounts_gemma_logfile}* | wc -l`
	if [ "${nsuccess}" = "${genetic_chunks}" ]; then
		echo "07e-gwas_cellcounts_mvlmm.sh completed successfully on all chromosomes."
	else
		echo "Problem: 07e-gwas_cellcounts_mvlmm.sh only complete for ${nsuccess} of ${genetic_chunks} batches"
		exit 1
	fi

}


sections=("01" "02" "03" "04" "05" "06" "07")
checkFirstArg "$1" "${sections[@]}"

actions=("check" "upload")
checkSecondArg "$2" "${actions[@]}"



echo "Checking log files for $1"
eval "check_logs_$1"


if [ "$2" = "upload" ]
then
	echo "Uploading"
fi




