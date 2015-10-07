#!/bin/bash

set -e
source config

R --no-save --args ${methylation_adjusted_pcs} ${methylation_adjusted_pcs_sq} < resources/methylation/methylation_matrixeqtl_format.R

