#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_08a_logfile})
print_version


if [ "${EWAS_phenotypes}" = "NULL" ]
then
	echo "No phenotypes have been specified."
else
	echo ""
	echo "Performing Height EWAS on all"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${methylation_adjusted}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		150 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_allindiv.RData \
		${section_08_dir}/qqplot_Height_allindiv \
		${section_08_dir}/Height_allindiv

	echo ""
	echo "Performing Height EWAS on children"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${methylation_adjusted}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		18 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_children.RData \
		${section_08_dir}/qqplot_Height_children
		${section_08_dir}/Height_children

	echo ""
	echo "Performing Height EWAS on adults"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${methylation_adjusted}.RData \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		18 \
		150 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_adults.RData \
		${section_08_dir}/qqplot_Height_adults
        ${section_08_dir}/Height_adults
fi

echo "Successfully performed height EWAS"
