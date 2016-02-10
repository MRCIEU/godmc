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
	echo "Performing BMI EWAS on all"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		0 \
		150 \
		20 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_BMI_allindiv.RData \
		${section_08_dir}/qqplot_BMI_allindiv

	echo ""
	echo "Performing BMI EWAS on children"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		0 \
		18 \
		20 \
		${nongenetic_meth_pcs} \
        ${section_08_dir}/results_BMI_children.RData \
		${section_08_dir}/qqplot_BMI_children

	echo ""
	echo "Performing BMI EWAS on adults"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		BMI \
		18 \
		150 \
		20 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_BMI_adults.RData \
		${section_08_dir}/qqplot_BMI_adults
fi

echo "Successfully performed BMI EWAS"
