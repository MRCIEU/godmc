#!/bin/bash

set -e
source config

R --no-save --args ${methylation_adjusted} ${meth_chunks} < ${home_directory}/resources/methylation/aggregate_chunks.R