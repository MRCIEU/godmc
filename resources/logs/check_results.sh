#!/usr/bin/env bash


check_results_01 () {

	if [ -f "${cohort_descriptives}" ]; then
		echo "Cohort descriptives file present"
	else
		echo "Cohort descriptives file absent. Please re-run."
		exit 1
	fi

	if [ -f "${methylation_summary}" ]; then
		echo "Methylation summary file present"
	else
		echo "Methylation summary file absent. Please re-run."
		exit 1
	fi

}


check_results_02 () {

	if [ -f "${allele_ref}" ]; then
		echo "Allele reference file present"
	else
		echo "Problem: Allele reference file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.frq.gz" ]; then
		echo "Allele frequency file present"
	else
		echo "Problem: Allele frequency file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.hwe.gz" ]; then
		echo "HWE file present"
	else
		echo "Problem: HWE file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.info.gz" ]; then
		echo "Imputation quality file present"
	else
		echo "Problem: Imputation quality file is absent"
		exit 1
	fi

}

check_results_03 () {

	echo "No results expected"
}

check_results_04 () {

	if [ -f "${cellcounts_plot}" ]; then
		echo "Cellcounts plot present"
	elif [ -f "${cellcounts_gwa}" ]; then
		echo "Cellcounts plot file not present. Please check 04a."
		exit 1
	else
		echo "Cellcounts GWAS not being performed. GWAS on cell counts will not be performed."
	fi
	if [ -f "${smoking_pred_plot}" ]; then
		echo "Smoking prediction plot present"
	else
		echo "Smoking prediction plot file not present"
		exit 1
	fi
	if [ -f "${age_pred_plot}" ]; then
		echo "Age prediction plot present"
	else
		echo "WARNING: Age prediction plot file not present. GWAS on age acceleration will not be performed."
	fi
	
	if [ -f "${section_04_dir}/positive_control_pcadjusted_${positive_control_cpg}.qassoc.gz" ]; then
		echo "PC adjusted positive control results present"
	else
		echo "PC adjusted positive control results file not present"
		exit 1
	fi
	if [ -f "${section_04_dir}/positive_control_pcadjusted_${positive_control_cpg}_manhattan.png" ]; then
		echo "PC adjusted positive control Manhattan plot present"
	else
		echo "PC adjusted positive control Manhattan plot file not present"
		exit 1
	fi
	
	if [ -f "${section_04_dir}/positive_control_pcadjusted_${positive_control_cpg}_qqplot.png" ]; then
		echo "PC adjusted positive control QQ plot present"
	else
		echo "PC adjusted positive control QQ plot file not present"
		exit 1
	fi
	if [ -f "${section_04_dir}/positive_control_pcunadjusted_${positive_control_cpg}.qassoc.gz" ]; then
		echo "PC unadjusted positive control results present"
	else
		echo "PC unadjusted positive control results file not present"
		exit 1
	fi
	if [ -f "${section_04_dir}/positive_control_pcunadjusted_${positive_control_cpg}_manhattan.png" ]; then
		echo "PC unadjusted positive control Manhattan plot present"
	else
		echo "PC unadjusted positive control Manhattan plot file not present"
		exit 1
	fi
	if [ -f "${section_04_dir}/positive_control_pcunadjusted_${positive_control_cpg}_qqplot.png" ]; then
		echo "PC unadjusted positive control QQ plot present"
	else
		echo "PC unadjusted positive control QQ plot file not present"
		exit 1
	fi

}

check_results_05 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`ls -l ${section_05_dir}/res* | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All mQTL results present"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} mQTL result files are present"
		exit 1
	fi

}

check_results_06 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`ls -l ${section_06_dir}/res* | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All variance-mQTL results present"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} variance-mQTL result files are present"
		exit 1
	fi

}


check_results_07 () {

	nbatch=`ls -l ${tabcnv}.tab.* | wc -l`
	nsuccess=`ls -l ${section_07_dir}/res* | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All CNV-mQTL results present"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} CNV-mQTL result files are present"
	fi

}



check_results_08 () {

	if [ "${phenotypes}" = "NULL" ]
	then
		echo "No phenotypes"
		exit 1
	fi

	if grep -q "Height" ${phenotype_list}; then
		if [ -f "${section_08_dir}/results_Height_allindiv.RData" ]; then
			echo "Height EWAS results present"
		else
			echo "Problem: Height EWAS results absent"
			exit 1
		fi
	fi

	if grep -q "BMI" ${phenotype_list}; then
		if [ -f "${section_08_dir}/results_BMI_allindiv.RData" ]; then
			echo "BMI EWAS results present"
		else
			echo "Problem: BMI EWAS results absent"
			exit 1
		fi
	fi
}


check_results_09 () {

	if grep -l -q "loco" ${section_09_dir}/*; then
		echo "AAR GWAS results present"
	else
		echo "Problem: AAR GWAS results absent"
		exit 1
	fi

}

check_results_10 () {

	if [ -f "${section_10_dir}/gwas_list.txt" ]; then
		echo "GWAS list annotation present"
	else
		echo "Problem: GWAS list is absent"
		exit 1
	fi

	IFS=$'\r\n' GLOBIGNORE='*' command eval  'pheno=($(cat ${section_10_dir}/gwas_list.txt))'

	for phen in ${pheno[@]}
	do
		if [ -f "${section_10_dir}/${phen}.loco.mlma.gz" ]; then
			echo "GWAS for ${phen} complete"
		else
			echo "problem: GWAS did not complete for ${phen}."
			exit 1
		fi
	done

}

check_results_11 () {

	if grep -l -q "loco" ${section_11_dir}/*; then
		echo "Cellcount entropy GWAS results present"
	else
		echo "Problem: Cellcount entropy GWAS results absent"
		exit 1
	fi

}


check_results_12 () {

	if [ -f "${section_12_dir}/cellcounts_columns.txt" ]; then
		echo "Cellcount annotation present"
	else
		echo "Problem: Cellcount annotation is absent"
		exit 1
	fi

	ncellcounts=`wc -l ${section_12_dir}/cellcounts_columns.txt | awk '{ print $1 }'`
	for i in $(seq 1 $ncellcounts);
	do
		if [ -f "${section_12_dir}/cellcount_${i}.loco.mlma.gz" ]; then
			echo "GWAS for cell type ${i} of ${ncellcounts} complete"
		else
			echo "problem: GWAS did not complete for cell type ${i}"
			exit 1
		fi
	done

}


check_results_13 () {

	nsuccess=`grep -l "assoc" ${section_13_dir}/* | wc -l`
	if [ "${nsuccess}" = "${genetic_chunks}" ]; then
		echo "Cellcounts results present for all chromosomes"
	else
		echo "Problem: Only ${nsuccess} of ${genetic_chunks} batches result files are present for the MVLMM cellcount GWAS"
		exit 1
	fi

}


check_results_14 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`ls -l ${section_14_dir}/res* | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All mQTL results present for X analysis on females"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} mQTL result files are present"
		exit 1
	fi

}

check_results_15 () {

	nbatch=`ls -l ${tabfile}.tab.* | wc -l`
	nsuccess=`ls -l ${section_15_dir}/res* | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All mQTL results present for X analysis on males"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} mQTL result files are present"
		exit 1
	fi

}

check_results_16 () {

	nbatch=`ls -l ${methylation_processed_dir}/methylation.subset*ge1.2.txt | wc -l`
	nsuccess=`ls -l ${section_16_dir}/gcta*.txt | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All GCTA results present"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} GCTA mQTL result files are present"
		exit 1
	fi

}check_results_17 () {

	nbatch="1"
	nsuccess=`ls -l ${section_17_dir}/plink*.txt | wc -l`
	if [ "${nsuccess}" = "${nbatch}" ]; then
		echo "All PLINK results present"
	else
		echo "Problem: Only ${nsuccess} of ${nbatch} PLINK mQTL result files are present"
		exit 1
	fi

}

