#!/bin/bash

# create batch scripts
# e.g to run:

# $ ./scriptname 1

set -e
source config

batch_number=${1}
re='^[0-9]+$'
if ! [[ $batch_number =~ $re ]] ; then
   echo "error: Batch variable is not a number"
   exit 1
fi

geno="${bfile}.${batch_number}.tab"
phen="${}"
threshold=${soft_threshold}
out="${}.${batch_number}.RData"

R --no-save --args ${geno} ${phen} ${cov} ${threshold} ${out} < ${matrixeqtl_run_r}


