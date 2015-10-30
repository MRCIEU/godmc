#!/bin/bash

set -e
source config
exec &> >(tee ${methylation_adjustment1_logfile}_aggregation)

echo "Aggregating methylation_adjustment1 chunks"
Rscript ${home_directory}/resources/methylation/aggregate_chunks.R \
	${methylation_adjusted} \
	${meth_chunks}