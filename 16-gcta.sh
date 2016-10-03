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

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

#cpgs2=`printf '${cpgs}\n%.0s' {1..$nocohorts}`

echo $batch_number
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

        probes=`grep -w $probe ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt | wc -l`

        if [ "$probes" -gt "0" ]   
        then
        
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $2}' > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
    
    #check whether fgrep is grepping one probe only
    #zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $3}' |sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.probe
    
    awk -F":" '{print $1}' < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps | sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

    Counter=`expr $Counter + 1`
    echo "$Counter/$noprobes probes"



            while read -r line
            do
                chr="$line"
                echo $chr
      
                grep -w $chr ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
                nosnps=`grep -f ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr ${bfile}.bim | wc -l`
                echo "no snps = $nosnps"

                if [ "$nosnps" -gt "0" ]
                then
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

                cat ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma | sed 's/^/'$probe'\t/'| perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
                tail -n +2 ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp >>${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt
        
                rm ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
                rm ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma

                fi
            
            done < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

            
        fi
        
done < "$filename"

echo "Successfully completed ${batch_number}"



