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
exec &> >(tee ${section_16_logfile}${batch_number})
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

noprobes=`cat $filename |wc -l`
echo $noprobes

echo "CpG" "CHR" "SNP" "BP" "NMISS" "BETA" "SE" "R2" "T" "P" |perl -pe 's/ /\t/g' >${lmm_res_dir}/plink.${i}\_${j}.ge${no}.txt

Counter=0

while read -r line
do
    probe="$line"
    echo $probe
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz |fgrep $probe |awk '{print $2}' > ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
  
    Counter=`expr $Counter + 1`
    echo "$Counter/$noprobes"

    nosnps=`grep -f ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps ${bfile}.bim |wc -l`
    echo $nosnps    
        if [ "$nosnps" -gt "0" ]
   
        then
     
        ${plink} \
            --bfile ${bfile} \
            --extract ${methylation_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps \
            --pheno ${methylation_processed_dir}/methylation.subset.${i}\_${j}.ge${no}.txt \
            --mpheno $Counter \
            --out ${lmm_res_dir}/plink.${i}\_${probe}.ge${no} \
            --assoc \
            --threads ${nthreads}



        sed 's/^/'$probe'/' <${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
        

        tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${lmm_res_dir}/plink.${i}\_${j}.ge${no}.txt
        
        rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
        rm ${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.qassoc
        rm ${lmm_res_dir}/plink.${i}\_${probe}.ge${no}.log
        
        fi
done < "$filename"

#merge effect alleles and AF to output files

echo "merge effect alleles and AF to output files"
cd ${lmm_res_dir}
sed 's/  \+/\t/g' < plink.${i}\_${j}.ge${no}.txt | sed 's/^ //g'  >plink.${i}\_${j}.ge${no}.txt.tmp
mv plink.${i}\_${j}.ge${no}.txt.tmp plink.${i}\_${j}.ge${no}.txt
cp plink.${i}\_${j}.ge${no}.txt plink.${i}.${PBS_ARRAYID}.ge${no}.txt
perl ${home_directory}/resources/phase2/join_file.pl -i "${lmm_res_dir}/plink.${i}.${PBS_ARRAYID}.ge${no}.txt,TAB,2 data.frq.tmp,TAB,1" -o ${lmm_res_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt -a 1
gzip ${lmm_res_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt

#rm plink.${i}.${PBS_ARRAYID}.ge${no}.txt
date

echo "Successfully completed"

