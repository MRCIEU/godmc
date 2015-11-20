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

	echo "No results expected"
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
		if [ -f "${section_08_dir}/results.Height.RData" ]; then
			echo "Height EWAS results present"
		else
			echo "Problem: Height EWAS results absent"
			exit 1
		fi
	fi

	if grep -q "BMI" ${phenotype_list}; then
		if [ -f "${section_08_dir}/results.BMI.RData" ]; then
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

	if grep -l -q "loco" ${section_10_dir}/*; then
		echo "Smoking GWAS results present"
	else
		echo "Problem: Smoking GWAS results absent"
		exit 1
	fi

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
