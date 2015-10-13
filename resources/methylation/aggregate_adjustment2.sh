#!/bin/bash

set -e
source config

R --no-save --args ${methylation_rt_cc} ${meth_chunks} < ${home_directory}/resources/methylation/aggregate_chunks.R