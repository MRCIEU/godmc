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
exec &> >(tee ${section_16b_logfile}${batch_number})
print_version

mydir="${methylation_processed_dir}"
cd $mydir

myout="${section_16_dir}"
echo $myout

i="1e-05"

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01" "cg02" "cg03" "cg04" "cg05" "cg06" "cg07" "cg08" "cg09" "cg10" "cg11" "cg12" "cg13" "cg14" "cg15" "cg16" "cg17" "cg18" "cg19" "cg20" "cg21" "cg22" "cg23" "cg24" "cg25" "cg26" "cg27" "_ch")
#cpgs2=`printf '${cpgs}\n%.0s' {1..$nocohorts}`

batch_number_zero=$((batch_number - 1))
j=${cpgs[${batch_number_zero}]}
echo $j

no="1.2"

filename=${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
echo $filename

noprobes=`cat $filename | wc -l`
echo $noprobes

echo "CpG" "CHR" "SNP" "BP" "EA" "NEA" "EAF" "BETA" "SE" "P" | perl -pe 's/ /\t/g' > ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt

Counter=0

while read -r line
do
    probe="$line"
    echo $probe
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $2}' > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
    
    #check fgrep
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $3}' |sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.probe
    
    awk -F":" '{print $1}' < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps | sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

    Counter=`expr $Counter + 1`
    echo "$Counter/$noprobes probes"

    nosnps=`grep -f ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps ${bfile}.bim | wc -l`
    echo "no snps= $nosnps"
        if [ "$nosnps" -gt "0" ]   
        then
            while read -r line
            do
                chr="$line"
                echo $chr
      
                grep -w $chr ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
                
                chrno=`echo $chr |sed 's/chr//g'`     

                ${plink} \
                    --bfile ${bfile} \
                    --extract ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --make-bed \
                    --out ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr

                if [ "$chrno" -lt "23" ]
                then
                ${gcta} \
                    --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --mlma \
                    --grm ${grmfile_all}_minus_chr${chrno} \
                    --pheno ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt \
                    --mpheno $Counter \
                    --qcovar ${covariates_combined}.gcta.numeric \
                    --covar ${covariates_combined}.gcta.factor \
                    --out ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno} \
                    --thread-num ${nthreads}
                fi
                    
                if [ "$chrno" -eq "23" ]
                then
                ${gcta} \
                    --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --mlma \
                    --grm ${grmfile_all} \
                    --pheno ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt \
                    --mpheno $Counter \
                    --qcovar ${covariates_combined}.gcta.numeric \
                    --covar ${covariates_combined}.gcta.factor \
                    --out ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno} \
                    --thread-num ${nthreads}
                fi        

            
            done < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

        cat ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma | sed 's/^/'$probe'\t/'| perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
        tail -n +2 ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp >>${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt
        
        rm ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
        rm ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma
        
        fi
        
done < "$filename"

echo "Successfully completed ${batch_number}"



