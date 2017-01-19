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

outfile="${section_16_dir}/batch_${batch_number}"

if [ -a "${outfile}" ]; then
    echo "Warning! Partial output file already exists:"
    echo "${outfile}"
    echo "Deleting old file before proceeding"
    rm ${outfile}
fi

if [ -a "${outfile}.gz" ]; then
    echo "Warning! Complete outfile already exists:"
    echo "${outfile}.gz"
    echo "Please delete this file to re-run the analysis"
    exit 1
fi

touch ${outfile}


i=0
nprobe=`cat ${phase2_assoclist}/assoclist_${batch_number}.probes | wc -l`

zcat ${phase2_assoclist}/assoclist_${batch_number}.gz | cut -d " " -f 2 | sort -u > ${phase2_assoclist}/assoclist_${batch_number}.snps

${plink} --noweb \
    --bfile ${bfile}_phase2 \
    --extract ${phase2_assoclist}/assoclist_${batch_number}.snps \
    --recode A \
    --out ${phase2_scratch}/geno_${batch_number}

gzip -f ${phase2_scratch}/geno_${batch_number}.raw

for probe in `cat ${phase2_assoclist}/assoclist_${batch_number}.probes`
do
    i=$((i + 1))

    echo ""
    echo ""
    echo "${i} of ${nprobe}"

    # Get SNP list from assoc list

    zgrep -w ${probe} ${phase2_assoclist}/assoclist_${batch_number}.gz | cat > ${phase2_scratch}/${probe}.list
    awk '{ if($4 == "c") { print $2 }}' ${phase2_scratch}/${probe}.list > ${phase2_scratch}/${probe}.cis_all
    awk '{ if($4 == "t") { print $2 }}' ${phase2_scratch}/${probe}.list > ${phase2_scratch}/${probe}.trans_all

    # Find SNPs present in current study

    zfgrep -wf ${phase2_scratch}/${probe}.cis_all ${section_16_dir}/snplist.txt.gz | cat > ${phase2_scratch}/${probe}.cis


    zfgrep -wf ${phase2_scratch}/${probe}.trans_all ${section_16_dir}/snplist.txt.gz | cat > ${phase2_scratch}/${probe}.trans


    ncis=`cat ${phase2_scratch}/${probe}.cis | wc -l`
    ntrans=`cat ${phase2_scratch}/${probe}.trans | wc -l`

    echo "  ${probe}: ${ncis} cis and ${ntrans} trans"

    # If there are any cis SNPs...

    if [ "$ncis" -gt "0" ]; then

        # Do GWAS

        ${plink} --noweb \
            --bfile ${bfile}_phase2 \
            --extract ${phase2_scratch}/${probe}.cis \
            --assoc \
            --allow-no-sex \
            --pheno ${phase2_betas}${batch_number} \
            --pheno-name ${probe} \
            --out ${phase2_scratch}/${probe}_cis 

        # Sort this file out
        # snp, nmiss, beta, se, p
        sed 1d ${phase2_scratch}/${probe}_cis.qassoc | awk -v probe="$probe" '{ print probe, $2, $4, $5, $6, $9, "c" }' > ${phase2_scratch}/${probe}_cis.txt
        cat ${phase2_scratch}/${probe}_cis.txt >> ${outfile}

    fi


    # If there are cis and trans SNPs, get the best cis SNP and make covariate

    if [ "$ncis" -gt "0" ] && [ "$ntrans" -gt "0" ]; then

        # Get best cis SNP and create covariate file

        sort -gk 6 ${phase2_scratch}/${probe}_cis.txt | head -n 1 | awk '{ print $2, $6 }' > ${phase2_scratch}/${probe}.bestcis
        echo "Best cis SNP:"
        cat ${phase2_scratch}/${probe}.bestcis

        bestcis=`awk '{ print $1 }' ${phase2_scratch}/${probe}.bestcis`
        ${plink} --noweb \
            --bfile ${bfile}_phase2 \
            --snp ${bestcis} \
            --recode A \
            --out ${phase2_scratch}/${probe}
        sed 1d ${phase2_scratch}/${probe}.raw | awk '{ print $1, $2, $7 }' > ${phase2_scratch}/${probe}.cov


        # Run linear regression on trans snps

        ${plink} --noweb \
            --bfile ${bfile}_phase2 \
            --extract ${phase2_scratch}/${probe}.trans \
            --linear \
            --allow-no-sex \
            --pheno ${phase2_betas}${batch_number} \
            --pheno-name ${probe} \
            --covar ${phase2_scratch}/${probe}.cov \
            --out ${phase2_scratch}/${probe}_trans


        # Sort this file out

        grep -w "ADD" ${phase2_scratch}/${probe}_trans.assoc.linear | awk -v probe="$probe" '{ if ($7 != "NA") { print probe, $2, $6, $7, $7/$8, $9, "t" }}' > ${phase2_scratch}/${probe}_trans.txt
        cat ${phase2_scratch}/${probe}_trans.txt >> ${outfile}

    fi


    # If there are only trans SNPs...

    if [ "$ncis" -eq "0" ] && [ "$ntrans" -gt "0" ]; then

        ${plink} --noweb \
            --bfile ${bfile}_phase2 \
            --extract ${phase2_scratch}/${probe}.trans \
            --assoc \
            --allow-no-sex \
            --pheno ${phase2_betas}${batch_number} \
            --pheno-name ${probe} \
            --out ${phase2_scratch}/${probe}_trans


        # Sort this file out
        sed 1d ${phase2_scratch}/${probe}_trans.qassoc | awk -v probe="$probe" '{ print probe, $2, $4, $5, $6, $9, "t" }' > ${phase2_scratch}/${probe}_trans.txt
        cat ${phase2_scratch}/${probe}_trans.txt >> ${outfile}

    fi

    rm -f ${phase2_scratch}/${probe}*

done


# Merge full results with SNP allele info etc

echo ""
echo ""
echo "Converting results to GWAMA format"

Rscript resources/phase2/create_gwama.R \
    ${outfile} \
    ${section_16_dir}/data.frq.gz

# Gzip
gzip ${outfile}


echo "Successfully completed"

