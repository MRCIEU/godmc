#!/bin/bash

set -e
source config
exec &> >(tee ${section_04b_logfile}_aggregation)

echo "Aggregating methylation_adjustment1 chunks"
Rscript ${home_directory}/resources/methylation/aggregate_chunks.R \
	${methylation_adjusted} \
	${meth_chunks}

Rscript ${home_directory}/resources/methylation/aggregate_chunk_classes.R \
	${methylation_adjusted} \
	${meth_chunks} \
	${section_04_dir}/adjustment1_classes.RData

rm ${methylation_adjusted}.*.RData