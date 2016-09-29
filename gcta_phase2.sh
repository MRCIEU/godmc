#!/bin/bash

#PBS -N plink
#PBS -o /panfs/panasas01/shared-godmc/job_report/plink.o
#PBS -e /panfs/panasas01/shared-godmc/job_report/plink.e
#PBS -l walltime=100:00:00
#PBS -t 1-74
# PBS -t 3
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

cd /panfs/panasas01/sscm/epzjlm/repo/godmc
source ./config


mydir="${methylation_processed_dir}"
cd $mydir

myout="${section_16_dir}"
echo $myout

i="1e-05"

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01" "cg02" "cg03" "cg04" "cg05" "cg06" "cg07" "cg08" "cg09" "cg10" "cg11" "cg12" "cg13" "cg14" "cg15" "cg16" "cg17" "cg18" "cg19" "cg20" "cg21" "cg22" "cg23" "cg24" "cg25" "cg26" "cg27" "_ch")
#cpgs2=`printf '${cpgs}\n%.0s' {1..$nocohorts}`

noa=$((PBS_ARRAYID - 1))

no=$((noa / ${#cpgs[@]}))
#no=$(($PBS_ARRAYID / ${#cpgs[@]}))

nob=$((${#cpgs[@]} * $no))
noc=$(($PBS_ARRAYID-$nob))
j=${cpgs[$noc-1]}
echo $j

no="1.2"

filename=${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
echo $filename

noprobes=`cat $filename |wc -l`
echo $noprobes

echo "CpG" "CHR" "SNP" "BP" "NMISS" "BETA" "SE" "R2" "T" "P" |perl -pe 's/ /\t/g' >${lmm_res_dir}/gcta.${i}\_${j}.ge${no}.txt

Counter=0

while read -r line
do
    probe="$line"
    echo $probe
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz |fgrep $probe |awk '{print $2}' > ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
    awk -F":" '{print $1}' <${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps | sort -u > ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

    Counter=`expr $Counter + 1`
    echo "$Counter/$noprobes probes"

    nosnps=`grep -f ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps ${bfile}.bim |wc -l`
    echo "no snps= $nosnps"    
        if [ "$nosnps" -gt "0" ]
   
        then
            
            while read -r line
            do
            chr="$line"
            echo $chr
  
            grep $chr ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps > ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
            
            chrno=`echo $chr |sed 's/chr//g'`     

            ${plink} \
            --bfile ${bfile} \
            --extract ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
            --make-bed \
            --out ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr


            ${gcta} \
            --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
            --mlma \
            --grm ${section_16_dir}/grm.minuschr${chrno} \
            --pheno ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt \
            --mpheno $Counter \
            --qcovar ${covariates_combined}.gcta.numeric \
            --covar ${covariates_combined}.gcta.factor \
            --chr $chrno \
            --out ${lmm_res_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.txt \
            --thread-num ${nthreads}
            
            done < ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs


        fi

        #sed 's/^/'$probe'/' <${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
        

        #tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${lmm_res_dir}/plink.${i}\_${j}.ge${no}.txt
        
        #rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
        #rm ${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.qassoc
        #rm ${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.log
        
        
done < "$filename"





