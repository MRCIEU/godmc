#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_13a_logfile})
print_version

if [ ! -f ${cellcounts_gwa} ]
then
	echo "No multivariate cellcounts GWA will be performed; please note we only run GWAS on Houseman estimates"
	exit 0
fi

if [ -f "${cellcounts_gwa}" ]
then
echo "Preparing inputfiles for GEMMA"

awk '{print $1,$2}' < ${bfile}.fam >$bfile.indiv

${gcta} \
		--grm ${grmfile_all} \
		--keep $bfile.indiv \
		--make-grm-bin \
		--out ${grmfile_all}.filtered \
		--thread-num ${nthreads}

Rscript resources/genetics/gemma_files.R \
	${grmfile_all}.filtered \
	${cellcounts_tf} \
	${section_12_dir}/cellcounts_columns.txt
fi

echo "Successfully formatted data for GEMMA"
