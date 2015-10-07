#!/bin/bash

set -e
source config

R --no-save --args ${methylation_adjusted_pcs}.RData ${normalised_phenotypes} ${ewas_results} < resources/methylation/ewas.R
