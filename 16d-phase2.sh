#!/bin/bash

set -e
source ./config

batch_number=${1}
re='^[0-9]+$'
if ! [[ $batch_number =~ $re ]] ; then
	echo "error: Batch variable is not a number"
	echo "Usage: ${0} [batch number]"
	exit 1
fi
# exec &> >(tee ${section_05_logfile}${batch_number})
print_version

n_genetic_batch=`ls -l ${tabfile}.tab.* | wc -l`
echo "Genetic data is split into ${n_genetic_batch} chunks"

if [ ! "${n_genetic_batch}" = "${genetic_chunks}" ]
then
	echo "Problem: Genetic data has been split into ${n_genetic_batch}, but the number of batches specified in the config file is ${genetic_chunks}"
	echo "Please either change the 'genetic_chunks' variable in the config file to ${n_genetic_batch} or re-run script 02b"
	exit
fi


geno="${tabfile}.tab.${batch_number}"
phen="${methylation_adjusted_pcs}.txt"
cov="NULL"
threshold="0.5"
out="temp.RData"

head -n 5001 ${geno} > tempgeno.tab
head -n 5001 ${phen} > tempphen.tab

geno="tempgeno.tab"
phen="tempphen.tab"


nbatch=(${tabfile}.tab.*)
nbatch=${#nbatch[@]}

echo "Performing meQTL analysis batch ${batch_number} of ${nbatch}"
Rscript resources/genetics/matrixeqtl.R ${geno} ${phen} ${cov} ${threshold} ${out}

echo "Successfully completed ${batch_number} of ${nbatch}"
