#!/bin/bash

set -e
source config

R --no-save --args ${methylation_rt_cc_poly} ${meth_chunks} < ${home_directory}/resources/methylation/aggregate_chunks.R