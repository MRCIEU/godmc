#!/bin/bash

set -e
source config

R --no-save --args ${betas} ${methylation_rt} ${nthreads} < resources/cellcounts/rank_transform.R

