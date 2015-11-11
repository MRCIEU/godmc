#!/bin/bash

set -e
source config
exec &> >(tee ${section_04d_logfile}_aggregation)

echo "Aggregating methylation_adjustment2 chunks"
Rscript ${home_directory}/resources/methylation/aggregate_chunks.R \
	${methylation_adjusted_pcs} \
	${meth_chunks}

rm ${methylation_adjusted_pcs}.*.RData