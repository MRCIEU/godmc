#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_01_logfile})
print_version


containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "./01-check_data [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg="all"
declare -a sections=('all' 'config' 'download' 'requirements' 'genetic' 'methylation' 'covariates' 'phenotypes' 'cnv' 'summary')

if [ -n "${1}" ]; then
	arg="${1}"
	containsElement ${1} ${sections[@]}
fi





section_message () {

	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ./01-check_data.sh $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""

}


if [ "$arg" = "config" ] || [ "$arg" = "all" ]
then

	section_message "config"

	if ! [[ "$study_name" =~ [^a-zA-Z0-9_\ ] ]] && ! [ "$study_name" = "" ]
	then
		echo "The study name '${study_name}' will be used for this analysis. Change this in the config file if necessary."
		echo ""
	else
		echo "The study name '${study_name}' is invalid. Please use only alphanumeric or underscore characters, with no spaces or special characters etc."
		exit
	fi

fi


if [ "$arg" = "download" ] || [ "$arg" = "all" ]
then

	section_message "download"



	sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources <<EOF
get 1kg_phase3_eur_allchrs_polymorphic.recoded.nodup.frq
get 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz.md5sum
EOF

	mv 1kg_phase3_eur_allchrs_polymorphic.recoded.nodup.frq.gz* ${home_directory}/resources/genetics

fi


if [ "$arg" = "requirements" ] || [ "$arg" = "all" ]
then

	section_message "requirements"


	Rscript resources/datacheck/requirements.R
fi

if [ "$arg" = "genetic" ] || [ "$arg" = "all" ]
then

	section_message "genetic"


	Rscript resources/datacheck/genetic_data.R \
		${bfile_raw}.bim \
		${bfile_raw}.fam \
		${quality_scores} \
		${control_snps} \
		${snpchrtxt} \
		${snpchrplot} \
		${genetic_descriptives} \
		${quality_scores_plot}

	# Check missingness, there should be a small percentage of missingness (--hard-call-threshold 0.499999)
	${plink} --bfile ${bfile_raw} --missing gz --out ${section_01_dir}/data

	nrow=`zcat ${section_01_dir}/data.imiss.gz | awk 'NR>1 && $6>0.02 {print $0}'  |wc -l`

	nrow_all=`zcat ${section_01_dir}/data.imiss.gz | awk 'NR>1 {print $0}' |wc -l`

	prop=$( printf '%.2f' $(echo "$nrow / $nrow_all" | bc -l) )
	echo "Proportion with >= 0.02 missingness: $prop"

	if [ "$prop" \> 0.01 ]
	then
		echo "Error: $prop of genotypes have more than 2% of missing values. Please don't use a genotype probability cut-off when generating best guess data."
		exit 1
	else
		echo "\n"
		echo "Best guess data appears to be correct"
	fi

fi

if [ "$arg" = "methylation" ] || [ "$arg" = "all" ]
then

	section_message "methylation"


	Rscript resources/datacheck/methylation_data.R \
		${betas} \
		${bfile_raw}.fam \
		${predicted_cellcounts} \
		${predicted_cellcounts_type} \
		${measured_cellcounts} \
		${meth_ids} \
		${methylation_descriptives} \
		${methylation_summary} \
		${intersect_ids} \
		${intersect_ids_plink}
fi

if [ "$arg" = "covariates" ] || [ "$arg" = "all" ]
then

	section_message "covariates"


	Rscript resources/datacheck/covariates.R \
		${covariates} \
		${bfile_raw}.fam \
		${meth_ids} \
		${ageplot} \
		${covariate_descriptives}
fi

if [ "$arg" = "phenotypes" ] || [ "$arg" = "all" ]
then

	section_message "phenotypes"


	Rscript resources/datacheck/phenotypes.R \
		${phenotypes} \
		${meth_ids} \
		${covariates} \
		${phenotype_descriptives} \
		${phenotype_list}
fi

if [ "$arg" = "cnv" ] || [ "$arg" = "all" ]
then

	section_message "cnv"


	Rscript resources/datacheck/cnv_data.R \
		${cnvs} \
		${meth_ids} \
		${covariates} \
		${cnv_descriptives}
fi

if [ "$arg" = "summary" ] || [ "$arg" = "all" ]
then

	section_message "summary"


	Rscript resources/datacheck/collect_descriptives.R \
		${genetic_descriptives} \
		${methylation_descriptives} \
		${covariate_descriptives} \
		${phenotype_descriptives} \
		${cnv_descriptives}	\
		${cohort_descriptives}
fi


if [ "$arg" = "all" ]
then
	echo ""
	echo ""
	echo "You successfully performed all data checks!"
fi
