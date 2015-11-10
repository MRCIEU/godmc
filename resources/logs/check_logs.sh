#!/usr/bin/env bash



check_logs_01 () {

	if grep -i -q "success" ${section_01_logfile}; then
		echo "01-check_data.sh completed successfully."
	else
		echo "Problem: 01-check_data.sh did not complete successfully"
		exit 1
	fi

}


check_logs_02 () {


	if grep -i -q "success" ${section_02a_logfile}; then
		echo "02a-snp_data.sh completed successfully."
	else
		echo "Problem: 02a-snp_data.sh did not complete successfully"
		echo "Please ensure previous steps have been performed and rerun 02a-snp_data.sh"
		exit 1
	fi


	if grep -i -q "success" ${section_02b_logfile}; then
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

	if grep -i -q "success" ${section_03a_logfile}; then
		echo "03a-phenotype_data.sh completed successfully."
	else
		echo "Problem: 03a-phenotype_data.sh did not complete successfully"
		exit 1
	fi

	if grep -i -q "success" ${section_03b_logfile}; then
		echo "03b-height_prediction.sh completed successfully."
	else
		echo "Problem: 03b-height_prediction.sh did not complete successfully"
		exit 1
	fi

}


check_logs_04 () {

	if grep -i -q "success" ${section_04a_logfile}; then
		echo "04a-methylation_variables.sh completed successfully."
	else
		echo "Problem: 04a-methylation_variables.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${section_04b_logfile}*; then
		echo "04b-methylation_adjustment1.sh completed successfully."
	else
		echo "Problem: 04b-methylation_adjustment1.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${section_04c_logfile}; then
		echo "04c-methylation_pcs.sh completed successfully."
	else
		echo "Problem: 04c-methylation_pcs.sh did not complete successfully"
		exit 1
	fi


	if grep -i -q "success" ${section_04d_logfile}*; then
		echo "04d-methylation_adjustment2.sh completed successfully."
	else
		echo "Problem: 04d-methylation_adjustment2.sh did not complete successfully"
		exit 1
	fi

	if grep -i -q "success" ${section_04e_logfile}*; then
		echo "04e-convert_methylation_format.sh completed successfully."
	else
		echo "Problem: 04e-convert_methylation_format.sh did not complete successfully"
		exit 1
	fi

}


check_logs_05 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`tail ${section_05_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "05-mqtl.sh completed successfully for all batches"
	else
		echo "Problem: 05-mqtl.sh only ${nsuccess} of ${nbatch} mQTL batches completed"
		exit 1
	fi

}

check_logs_06 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`tail ${section_06_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "06-vmqtl.sh completed successfully for all batches"
	else
		echo "Problem: 06-vmqtl.sh only ${nsuccess} of ${nbatch} variance-mQTL batches completed"
		exit 1
	fi

}

check_logs_07 () {

	nbatch=`ls -l ${tabcnv}.tab.* | wc -l`
	nsuccess=`tail ${section_07_logfile}* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "07-mcnv.sh completed successfully for all batches"
	else
		echo "Problem: 07-mcnv.sh only ${nsuccess} of ${nbatch} CNV-mQTL batches completed"
		exit 1
	fi

}


check_logs_08 () {

	if [ "${phenotypes}" = "NULL" ]
	then
		echo "No phenotypes specified"
		exit 1
	fi

	if grep -i -q "success" ${section_08_logfile}; then
		echo "08-ewas.sh completed successfully."
	else
		echo "Problem: 08-ewas.sh did not complete successfully"
		exit 1
	fi

}



check_logs_09 () {

	if grep -i -q "success" ${section_09_logfile}; then
		echo "09-gwas_aar.sh completed successfully."
	else
		echo "Problem: 09-gwas_aar.sh did not complete successfully"
		exit 1
	fi

}

check_logs_10 () {

	if grep -i -q "success" ${section_10_logfile}; then
		echo "10-gwas_smoking.sh completed successfully."
	else
		echo "Problem: 10-gwas_smoking.sh did not complete successfully"
		exit 1
	fi

}

check_logs_11 () {

	if grep -i -q "success" ${section_11_logfile}; then
		echo "11-gwas_cellcount_entropy.sh completed successfully."
	else
		echo "Problem: 11-gwas_cellcount_entropy.sh did not complete successfully"
		exit 1
	fi

}

check_logs_12 () {

	nbatch=`wc -l ${gwas_cellcounts_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
	nsuccess=`tail ${section_12_logfile}_* | grep -i "success" | wc -l`
	if [ "${nbatch}" = "${nsuccess}" ]; then
		echo "12-gwas_cellcounts.sh completed successfully for all cell types"
	else
		echo "problem: 12-gwas_cellcounts.sh only complete for ${nsuccess} of ${nbatch} cell types"
		exit 1
	fi

}

check_logs_13 () {

	nsuccess=`grep -i "success" ${section_13_logfile}* | wc -l`
	if [ "${nsuccess}" = "${genetic_chunks}" ]; then
		echo "13-gwas_cellcounts_mvlmm.sh completed successfully on all chromosomes."
	else
		echo "Problem: 13-gwas_cellcounts_mvlmm.sh only complete for ${nsuccess} of ${genetic_chunks} batches"
		exit 1
	fi

}