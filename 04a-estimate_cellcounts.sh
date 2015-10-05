#!/bin/bash

set -e
source config


# Estimate cell counts
if [ "${cellcounts_required}" = "yes" ]
then

	if [ ! "${provided_cellcounts}" = "NULL" ]
	then
		R --no-save --args ${betas} ${cellcounts} ${cellcount_reference} < resources/cellcounts/estimate_cellcounts.R
	elif [ -f "${provided_cellcounts}" ]
	then
		echo "Using the cellcounts provided in ${provided_cellcounts}."
	else
		echo "Error: The file ${provided_cellcounts} doesn't exist. You have specified that cell counts are required. Please set 'provided_cellcounts' to NULL if you want them to be estimated now, or specify a path to a file with the pre-specified cell counts."
	fi
elif [ "${cellcounts_required}" = "no" ]
then
	cellcounts="NULL"
else
	echo "'cellcounts_required' should be set to yes or no"
	exit 1
fi

