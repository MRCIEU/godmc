#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16b_logfile})
print_version

# Get the control CPG - in parameters file

# 1. Do GWAS on cg07959070 using just the cis markers
# 2. Keep the best one, create covariate
# 3. Check that it has p < 0.001
# 4. Create covariate, do GWAS and create manhattan plot, calculate lambda

echo "Running positive control"
echo "Testing cis-region of "

probe="cg07959070"
batch_number="202"
probe_chr="22"
probe_pos="50053871"
window_size="100000"
pval_threshold="0.001"



echo ""
echo ""
echo ""
echo "Running positive control - This involves testing the cis-region of ${probe}, extracting the top hit and performing GWAS with the cis-hit fitted as a covariate"
echo ""
echo "There should be a strong signal at the cis region, and no evidence for population stratification in the GWAS"
echo ""
echo ""
echo ""

# Check that control probe is present

pres=`zgrep -w ${probe} ${phase2_betas}${batch_number} | wc -l`

if [ "$pres" -eq "1" ]; then
	echo "Probe is available in the data"
else
	echo ""
	echo "WARNING: Probe is not present in the data."
	echo "WARNING: Please contact developers before continuing"
	# exit 1
fi

mkdir -p ${section_16_dir}/control


# Get SNP list from assoc list

zgrep -w ${probe} ${phase2_assoclist}/assoclist_${batch_number}.gz | cat > ${phase2_scratch}/${probe}.list



awk '{ if($4 == "c") { print $2 }}' ${phase2_scratch}/${probe}.list > ${phase2_scratch}/${probe}.cis_all


# Find SNPs present in current study

zfgrep -wf ${phase2_scratch}/${probe}.cis_all ${section_16_dir}/snplist.txt.gz | cat > ${phase2_scratch}/${probe}.cis


ncis=`cat ${phase2_scratch}/${probe}.cis | wc -l`

if [ "$ncis" -gt "0" ]; then
	echo "${probe}: ${ncis} cis SNPs"
else 
	echo "error: No cis SNPs identified. This is a problem, please contact developers."
	exit 1
fi

# Do GWAS on cis SNPs

${plink} --noweb \
	--bfile ${bfile}_phase2 \
	--extract ${phase2_scratch}/${probe}.cis \
	--assoc \
	--allow-no-sex \
	--pheno ${phase2_betas}${batch_number} \
	--pheno-name ${probe} \
	--out ${phase2_scratch}/${probe}_cis 

sed 1d ${phase2_scratch}/${probe}_cis.qassoc | awk -v probe="$probe" '{ print probe, $2, $4, $5, $6, $9, "c" }' | sort -gk 6 > ${section_16_dir}/control/${probe}_cis.txt


# Is the top hit p-value < threshold?

nsig=`awk -v thresh="${pval_threshold}" '{ if($6 < thresh) { print $0 }}' ${section_16_dir}/control/${probe}_cis.txt | wc -l`

head -n 1 ${section_16_dir}/control/${probe}_cis.txt | awk '{ print $2, $6 }' > ${phase2_scratch}/${probe}.bestcis
echo "Best cis SNP:"
cat ${phase2_scratch}/${probe}.bestcis

echo "${nsig} SNP(s) with p-value < ${pval_threshold}"

if [ "$nsig" -eq "0" ]; then
	echo ""
	echo ""
	echo "WARNING: Control probe has not found a strong association"
	echo "WARNING: Please contact developers"
fi




bestcis=`awk '{ print $1 }' ${phase2_scratch}/${probe}.bestcis`
${plink} --noweb \
	--bfile ${bfile}_phase2 \
	--snp ${bestcis} \
	--recode A \
	--out ${phase2_scratch}/${probe}
sed 1d ${phase2_scratch}/${probe}.raw | awk '{ print $1, $2, $7 }' > ${phase2_scratch}/${probe}.cov

echo "Running GWAS with cis-SNP fitted as covariate"

${plink} --noweb \
	--bfile ${bfile}_phase2 \
	--linear \
	--allow-no-sex \
	--pheno ${phase2_betas}${batch_number} \
	--pheno-name ${probe} \
	--covar ${phase2_scratch}/${probe}.cov \
	--out ${section_16_dir}/control/${probe}

rm -f ${section_16_dir}/control/${probe}.assoc.linear.gz

grep -w "ADD" ${section_16_dir}/control/${probe}.assoc.linear | gzip -c > ${section_16_dir}/control/${probe}.assoc.linear.gz

rm ${section_16_dir}/control/${probe}.assoc.linear


Rscript resources/genetics/plot_gwas.R \
	${section_16_dir}/control/${probe}.assoc.linear.gz \
	9 \
	1 \
	3 \
	FALSE \
	${probe_chr} \
	${probe_pos} \
	${window_size} \
	${pval_threshold} \
	${section_16_dir}/control/${probe}



echo "Successfully completed script 16b"

