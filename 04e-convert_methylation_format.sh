#!/bin/bash

set -e
source config
exec &> >(tee ${convert_methylation_format_logfile})

Rscript resources/methylation/methylation_matrixeqtl_format.R \
	${methylation_adjusted_pcs} \
	${methylation_adjusted_pcs_sq}

echo "Successfully converted methylation data"
