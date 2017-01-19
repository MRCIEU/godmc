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
exec &> >(tee ${section_16c_logfile}${batch_number})
print_version


# 1. Get batch
# 2. For each probe get cis and trans SNPs
# 3. Do GWAS on cis SNPs
# 4. Get best cis SNP and add to covariate file
# 5. Do GWAS on trans SNPs
# 6. Remove covariate file
# 7. Convert to gwama format and gzip


# Check if batch number is ok

nbatch=`ls -l ${phase2_betas}* | wc -l`
if [ "$batch_number" -gt "$nbatch" ] || [ "$batch_number" -lt "1" ]; then
    echo "error: Batch number must be between 1 and ${nbatch}"
    echo "Usage: ${0} [batch number]"
    exit 1
fi

echo "Performing section ${batch_number} of ${nbatch}"

outfile="${section_16_dir}/results_${batch_number}.gz"

if [ -a "${outfile}" ]; then
    echo "Warning! Output file already exists:"
    echo "${outfile}"
    echo "This will be overwritten"
fi

zcat ${phase2_assoclist}/assoclist_${batch_number}.gz | cut -d " " -f 2 | sort -u > ${phase2_scratch}/geno_${batch_number}.snplist

${plink} --noweb \
    --bfile ${bfile}_phase2 \
    --extract ${phase2_scratch}/geno_${batch_number}.snplist \
    --recode A \
    --out ${phase2_scratch}/geno_${batch_number}

gzip -f ${phase2_scratch}/geno_${batch_number}.raw
rm ${phase2_scratch}/geno_${batch_number}.snplist

Rscript resources/phase2/run_analysis.R \
    ${phase2_scratch}/geno_${batch_number}.raw.gz \
    ${phase2_betas}${batch_number} \
    ${phase2_assoclist}/assoclist_${batch_number}.gz \
    ${section_16_dir}/data.frq.gz \
    ${outfile}

rm ${phase2_scratch}/geno_${batch_number}.gz

echo "Successfully completed batch"
