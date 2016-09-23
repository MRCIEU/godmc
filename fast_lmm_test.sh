#!/bin/bash

#PBS -N clump
#PBS -o /panfs/panasas01/shared-godmc/job_report/clump
#PBS -e /panfs/panasas01/shared-godmc/job_report/clump
#PBS -l walltime=100:00:00
#PBS -t 1-47
# PBS -t 3
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi


mydir="/panfs/panasas01/shared-godmc/counts"

pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")
#pvals=("1e-13")

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

i="1e-05"
j="cg0000[0-9]"
no="1"
#cat $mydir/combined/cis.${i}\_${j}.allcohorts.txt $mydir/combined/trans.${i}\_${j}.allcohorts.txt > $mydir/combined/cis_trans.${i}\_${j}.allcohorts.txt
#perl -pe 's/_/ /g' $mydir/combined/cis_trans.${i}\_${j}.allcohorts.txt | awk '{print $3}' |sort -u > $mydir/combined/cis_trans.${i}\_${j}.allcohorts.probes

#/panfs/panasas01/sscm/epzjlm/repo/godmc/processed_data/methylation_data/cis_trans.1e-05_cg0000[0-9].ge1.allcohorts.probes
#/panfs/panasas01/sscm/epzjlm/repo/godmc/processed_data/methylation_data/cis_trans.1e-05\_cg0000[0-9].ge1.allcohorts.probes
#
#noqcov=`(head -n1 ${covariates_combined}.boltlmm |awk '{print gsub("QCOV","")}')`
#nocatcov=`(head -n1 ${covariates_combined}.boltlmm |awk '{print gsub("CATCOV","")}')`
#zcat ${hm3_snps_no_ld} >${hm3_snps_no_ld2}

filename=${methylation_processed_dir}/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
echo $filename


while read -r line
do
    probe="$line"
    echo $probe
fgrep $probe ${methylation_processed_dir}/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt |awk '{print $2}' > ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps

fastlmmc \
  -bfile ${bfile} \
  -extract ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps \
  -pheno ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt \
  -mpheno 1 \
  -out ${boltlmm_res_dir}/fastlmm.${i}\_${probe}.ge${no}.txt \
  -covar ${covariates_combined}.fastlmm \
  -sim ${grmfile_all} \
  -maxThreads 16

done < $filename





