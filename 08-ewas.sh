#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_08_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else
	echo ""
	echo "Performing Height EWAS on all"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		150 \
		${section_08_dir}/results_Height_all.RData \
		${section_08_dir}/qqplot_Height_all.png

	echo ""
	echo "Performing BMI EWAS on all"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		0 \
		150 \
		${section_08_dir}/results_BMI_all.RData \
		${section_08_dir}/qqplot_BMI_all.png

	echo ""
	echo "Performing Height EWAS on children"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		18 \
		${section_08_dir}/results_Height_children.RData \
		${section_08_dir}/qqplot_Height_children.png

	echo ""
	echo "Performing BMI EWAS on children"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		0 \
		18 \
		${section_08_dir}/results_BMI_children.RData \
		${section_08_dir}/qqplot_BMI_children.png

	echo ""
	echo "Performing Height EWAS on adults"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		18 \
		150 \
		${section_08_dir}/results_Height_adults.RData \
		${section_08_dir}/qqplot_Height_adults.png

	echo ""
	echo "Performing BMI EWAS on adults"
	echo ""
	Rscript resources/methylation/ewas.R \
		${methylation_adjusted_pcs}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		18 \
		150 \
		${section_08_dir}/results_BMI_adults.RData \
		${section_08_dir}/qqplot_BMI_adults.png
fi

echo "Successfully performed EWAS"
