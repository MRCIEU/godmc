#!/bin/bash

set -e
source config
exec &> >(tee ${gwas_cellcounts_logfile})

nval=`awk '{ print NR }' ${gwas_cellcounts_dir}/cellcounts_columns.txt | tr '\n' ' '`

${gemma} -bfile temp \
	-k ${grmfile_all}.gemma \
	-p ${cellcounts_tf}.gemma \
	-n 1 2 -lmm -o cellcounts_mvlmm


${gemma} -bfile temp \
	-k ${grmfile_all}.gemma \
	-p ~/sandpit/GEMMA/examplemouse_hs1940.pheno.txt \
	-n 1 2 -lmm -o cellcounts_mvlmm


${gwas_cellcounts_dir}/
echo "Successfully performed multivariate LMM on cellcounts"
