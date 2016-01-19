#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_01_logfile})
print_version

arg="all"

if [ -n "${1}" ]; then
	arg="${1}"
fi

if [ "$arg" = "download" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "download section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh download"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""



sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources <<EOF
get 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz
get 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz.md5sum
EOF
mv 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz* ${home_directory}/resources/genetics

fi


if [ "$arg" = "requirements" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "requirements section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh requirements"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


	Rscript resources/datacheck/requirements.R
fi

if [ "$arg" = "genetic" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "genetic section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh genetic"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


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

echo "-----------------------------------------------"
echo ""
echo "methylation section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh methylation"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


	Rscript resources/datacheck/methylation_data.R \
		${betas} \
		${bfile_raw}.fam \
		${provided_cellcounts} \
		${meth_ids} \
		${methylation_descriptives} \
		${methylation_summary} \
		${intersect_ids} \
		${intersect_ids_plink}
fi

if [ "$arg" = "covariates" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "covariates section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh covariates"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


	Rscript resources/datacheck/covariates.R \
		${covariates} \
		${bfile_raw}.fam \
		${meth_ids} \
		${ageplot} \
		${covariate_descriptives}
fi

if [ "$arg" = "phenotypes" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "phenotypes section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh phenotypes"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


	Rscript resources/datacheck/phenotypes.R \
		${phenotypes} \
		${meth_ids} \
		${covariates} \
		${phenotype_descriptives} \
		${phenotype_list}
fi

if [ "$arg" = "cnv" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "cnv section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh cnv"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


	Rscript resources/datacheck/cnv_data.R \
		${cnvs} \
		${meth_ids} \
		${covariates} \
		${cnv_descriptives}
fi

if [ "$arg" = "summary" ] || [ "$arg" = "all" ]
then

echo "-----------------------------------------------"
echo ""
echo "summary section"
echo ""
echo "to run this part on its own type:"
echo "$ ./01-check_data.sh summary"
echo ""
echo "-----------------------------------------------"
echo ""
echo ""


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
