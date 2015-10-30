#!/bin/bash

set -e
source config
exec &> >(tee ${methylation_adjustment2_logfile}_aggregation)

echo "Aggregating methylation_adjustment2 chunks"
Rscript ${home_directory}/resources/methylation/aggregate_chunks.R \
	${methylation_adjusted_pcs} \
	${meth_chunks}