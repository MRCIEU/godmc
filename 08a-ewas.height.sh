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
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		150 \
		20 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_allindiv.RData \
		${section_08_dir}/qqplot_Height_allindiv

	echo ""
	echo "Performing Height EWAS on children"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		0 \
		18 \
		20 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_children.RData \
		${section_08_dir}/qqplot_Height_children

	echo ""
	echo "Performing Height EWAS on adults"
	echo ""
	Rscript resources/methylation/ewas.meffil.R \
		${betas} \
		${ewastransformed} \
		${covariates_combined}.txt \
		Height \
		18 \
		150 \
		20 \
		${nongenetic_meth_pcs} \
		${section_08_dir}/results_Height_adults.RData \
		${section_08_dir}/qqplot_Height_adults

fi

echo "Successfully performed height EWAS"
