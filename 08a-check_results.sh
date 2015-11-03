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



if grep -i -q "success" ${snp_data_logfile}; then
	echo "02a-snp_data.sh completed successfully."
else
	echo "Problem: 02a-snp_data.sh did not complete successfully"
fi


if grep -i -q "success" ${convert_snp_format_logfile}; then
	echo "02b-convert_snp_format.sh completed successfully."
else
	echo "Problem: 02b-convert_snp_format.sh did not complete successfully"
fi


if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -i -q "success" ${phenotype_data_logfile}; then
		echo "03a-phenotype_data.sh completed successfully."
	else
		echo "Problem: 03a-phenotype_data.sh did not complete successfully"
	fi
fi


if [ ! "${phenotypes}" = "NULL" ]
then
	if grep -i -q "success" ${height_prediction_logfile}; then
		echo "03b-height_prediction.sh completed successfully."
	else
		echo "Problem: 03b-height_prediction.sh did not complete successfully"
	fi
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


nbatch=`ls -l ${tabfile}.tab.* | wc -l`
nsuccess=`tail ${mqtl_logfile}* | grep -i "success" | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "05a-mqtl.sh completed successfully for all batches"
else
	echo "Problem: 05a-mqtl.sh only ${nsuccess} of ${nbatch} mQTL batches completed"
fi


nbatch=`ls -l ${tabfile}.tab.* | wc -l`
nsuccess=`tail ${vmqtl_logfile}* | grep -i "success" | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "05b-vmqtl.sh completed successfully for all batches"
else
	echo "Problem: 05b-vmqtl.sh only ${nsuccess} of ${nbatch} variance-mQTL batches completed"
fi


nbatch=`ls -l ${tabcnv}.tab.* | wc -l`
nsuccess=`tail ${mcnv_logfile}* | grep -i "success" | wc -l`
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

nbatch=`wc -l ${gwas_cellcounts_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
nsuccess=`tail ${gwas_cellcounts_logfile}_* | grep -i "success" | wc -l`
if [ "${nbatch}" = "${nsuccess}" ]; then
	echo "07d-gwas_cellcounts.sh completed successfully for all cell types"
else
	echo "problem: 07d-gwas_cellcounts.sh only complete for ${nsuccess} of ${nbatch} cell types"
fi

nsuccess=`grep -i "success" ${gwas_cellcounts_gemma_logfile}* | wc -l`
if [ "${nsuccess}" = "${genetic_chunks}" ]; then
	echo "07e-gwas_cellcounts_mvlmm.sh completed successfully on all chromosomes."
else
	echo "Problem: 07e-gwas_cellcounts_mvlmm.sh only complete for ${nsuccess} of ${genetic_chunks} batches"
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

ncellcounts=`wc -l ${gwas_cellcounts_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
for i in $(seq 1 $ncellcounts);
do
	if [ -f "${gwas_cellcounts_dir}/cellcount_${i}.loco.mlma.gz" ]; then
		echo "GWAS for cell type ${i} of ${ncellcounts} complete"
	else
		echo "problem: GWAS did not complete for cell type ${i}"
	fi
done

nsuccess=`grep -l "assoc" ${gwas_cellcounts_dir}/* | wc -l`
if [ "${nsuccess}" = "${genetic_chunks}" ]; then
	echo "Cellcounts results present for all chromosomes"
else
	echo "Problem: Only ${nsuccess} of ${genetic_chunks} batches result files are present for the MVLMM cellcount GWAS"
fi

# mqtl
nbatch=`ls -l ${tabfile}.tab.* | wc -l`
nsuccess=`ls -l ${matrixeqtl_mqtl_dir}/res* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} mQTL result files are present"
fi

# SNP summaries
if [ -f "${allele_ref}" ]; then
	echo "Allele reference file present"
else
	echo "Problem: Allele reference file is absent"
fi

if [ -f "${matrixeqtl_mqtl_dir}/data.frq.gz" ]; then
	echo "Allele frequency file present"
else
	echo "Problem: Allele frequency file is absent"
fi

if [ -f "${matrixeqtl_mqtl_dir}/data.hwe.gz" ]; then
	echo "HWE file present"
else
	echo "Problem: HWE file is absent"
fi

if [ -f "${matrixeqtl_mqtl_dir}/data.info.gz" ]; then
	echo "Imputation quality file present"
else
	echo "Problem: Imputation quality file is absent"
fi


# vmqtl
nbatch=`ls -l ${tabfile}.tab.* | wc -l`
nsuccess=`ls -l ${matrixeqtl_vmqtl_dir}/res* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All variance-mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} variance-mQTL result files are present"
fi

# mcnv
nbatch=`ls -l ${tabcnv}.tab.* | wc -l`
nsuccess=`ls -l ${matrixeqtl_mcnv_dir}/res* | wc -l`
if [ "${nsuccess}" = "${nbatch}" ]; then
	echo "All CNV-mQTL results present"
else
	echo "Problem: Only ${nsuccess} of ${nbatch} CNV-mQTL result files are present"
fi

