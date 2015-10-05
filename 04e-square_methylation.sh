#!/bin/bash

set -e
source config

R --no-save --args ${methylation_rt_cc_poly} ${methylation_rt_cc_poly_sq} < resources/cellcounts/square_methdata.R

