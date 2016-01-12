#!/bin/bash

set -e
source config
exec &> >(tee ${section_01_logfile})
print_version


sftp ${sftp_username}@${sftp_address}:/GoDMC/resources <<EOF
get 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz
EOF
mv 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz ./resources/genetics

Rscript resources/datacheck/requirements.R

Rscript resources/datacheck/genetic_data.R \
	${bfile_raw}.bim \
	${bfile_raw}.fam \
	${quality_scores} \
	${control_snps} \
	${snpchrtxt} \
	${snpchrplot} \
	${genetic_descriptives} \
	${quality_scores_plot}

Rscript resources/datacheck/methylation_data.R \
	${betas} \
	${bfile_raw}.fam \
	${provided_cellcounts} \
	${meth_ids} \
	${methylation_descriptives} \
	${methylation_summary} \
	${intersect_ids} \
	${intersect_ids_plink}

Rscript resources/datacheck/covariates.R \
	${covariates} \
	${bfile_raw}.fam \
	${meth_ids} \
	${ageplot} \
	${covariate_descriptives}

Rscript resources/datacheck/phenotypes.R \
	${phenotypes} \
	${meth_ids} \
	${covariates} \
	${phenotype_descriptives} \
	${phenotype_list}

Rscript resources/datacheck/cnv_data.R \
	${cnvs} \
	${meth_ids} \
	${covariates} \
	${cnv_descriptives}

Rscript resources/datacheck/collect_descriptives.R \
	${genetic_descriptives} \
	${methylation_descriptives} \
	${covariate_descriptives} \
	${phenotype_descriptives} \
	${cnv_descriptives}	\
	${cohort_descriptives}

echo "You successfully performed the first data check"

# Check missingness, there should be a small percentage of missingness (--hard-call-threshold 0.499999)
${plink} --bfile ${bfile_raw} --missing gz --out ${section_01_dir}/data
 
nrow=`zcat ${section_01_dir}/data.imiss.gz | awk 'NR>1 && $6>0.01 {print $0}'  |wc -l`

nrow_all=`zcat ${section_01_dir}/data.imiss.gz | awk 'NR>1 {print $0}' |wc -l`

prop=$( printf '%.2f' $(echo "$nrow / $nrow_all" | bc -l) )
echo $prop


if [ "$prop" \> 0.01 ]
then
   echo "Error: $prop percent of genotypes have more than 1% of missing values. Please don't use a genotype probability cut-off."
   exit 1
 else
   echo "You successfully converted genotype data to bestguess format"
 fi
   echo "You successfully performed all data checks"

