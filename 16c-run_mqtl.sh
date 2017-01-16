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
exec &> >(tee ${section_16d_logfile}${batch_number})
print_version

mydir="${methylation_processed_dir}"
cd $mydir

myout="${section_16_dir}"
echo $myout

i="1e-05"

#subset 74
#cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

#subset 128
cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-2]" "cg01[3-5]" "cg01[6-8]" "cg01[9]" "cg02[0-2]" "cg02[3-5]" "cg02[6-8]" "cg02[9]" "cg03[0-2]" "cg03[3-5]" "cg03[6-8]" "cg03[9]" "cg04[0-2]" "cg04[3-5]" "cg04[6-8]" "cg04[9]" "cg05[0-2]" "cg05[3-5]" "cg05[6-8]" "cg05[9]" "cg06[0-2]" "cg06[3-5]" "cg06[6-8]" "cg06[9]" "cg07[0-2]" "cg07[3-5]" "cg07[6-8]" "cg07[9]" "cg08[0-2]" "cg08[3-5]" "cg08[6-8]" "cg08[9]" "cg09[0-2]" "cg09[3-5]" "cg09[6-8]" "cg09[9]" "cg10[0-2]" "cg10[3-5]" "cg10[6-8]" "cg10[9]" "cg11[0-2]" "cg11[3-5]" "cg11[6-8]" "cg11[9]" "cg12[0-2]" "cg12[3-5]" "cg12[6-8]" "cg12[9]" "cg13[0-2]" "cg13[3-5]" "cg13[6-8]" "cg13[9]" "cg14[0-2]" "cg14[3-5]" "cg14[6-8]" "cg14[9]" "cg15[0-2]" "cg15[3-5]" "cg15[6-8]" "cg15[9]" "cg16[0-2]" "cg16[3-5]" "cg16[6-8]" "cg16[9]" "cg17[0-2]" "cg17[3-5]" "cg17[6-8]" "cg17[9]" "cg18[0-2]" "cg18[3-5]" "cg18[6-8]" "cg18[9]" "cg19[0-2]" "cg19[3-5]" "cg19[6-8]" "cg19[9]" "cg20[0-2]" "cg20[3-5]" "cg20[6-8]" "cg20[9]" "cg21[0-2]" "cg21[3-5]" "cg21[6-8]" "cg21[9]" "cg22[0-2]" "cg22[3-5]" "cg22[6-8]" "cg22[9]" "cg23[0-2]" "cg23[3-5]" "cg23[6-8]" "cg23[9]" "cg24[0-2]" "cg24[3-5]" "cg24[6-8]" "cg24[9]" "cg25[0-2]" "cg25[3-5]" "cg25[6-8]" "cg25[9]" "cg26[0-2]" "cg26[3-5]" "cg26[6-8]" "cg26[9]" "cg27[0-2]" "cg27[3-5]" "cg27[6-8]" "cg27[9]" "_ch")


#cpgs2=`printf '${cpgs}\n%.0s' {1..$nocohorts}`

batch_number_zero=$((batch_number - 1))
j=${cpgs[${batch_number_zero}]}
echo $j

no="1.2"

filename=${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
echo $filename

noprobes=`cat $filename | wc -l`
echo $noprobes

echo "CpG" "CHR" "SNP" "BP" "NMISS" "BETA" "SE" "R2" "T" "P" |perl -pe 's/ /\t/g' >${section_16_dir}/plink.${i}\_${j}.ge${no}.txt

nprobes=`awk -F' ' '{print NF; exit}'< ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt`
nprobes=$((nprobes - 2))

echo "$nprobes/$noprobes probes found in methylation data"

Counter=0

while read -r line
do
    probe="$line"
    echo $probe
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz |fgrep $probe |awk '{print $2}' > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps

        probes=`grep -w $probe ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt | wc -l`

        if [ "$probes" -gt "0" ]   
        then
       
        Counter=`expr $Counter + 1`
        echo "$Counter/$noprobes"

            nosnps=`grep -f ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps ${bfile}.bim | wc -l`
            echo "no snps= $nosnps"
            
            #check whether right probes are selected from methylation matrix
            probecheck=`expr $Counter + 2`
            probecheck2=`head -n1 ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt | cut -f $probecheck -d" "` 

            echo "$probe is $probecheck2"
            
            if ! [[ "$probe" -eq "$probecheck2" ]] ; then
            echo "error: $probe doesn't match to methylation $probecheck2"
            exit 1
            fi
            
                if [ "$nosnps" -gt "0" ]
                then    

            ${plink} \
                --bfile ${bfile} \
                --extract ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps \
                --pheno ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt \
                --mpheno $Counter \
                --out ${section_16_dir}/plink.${i}\_${probe}.ge${no} \
                --allow-no-sex \
                --assoc \
                --threads ${nthreads}

            sed 's/^/'$probe'/' <${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
            tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${section_16_dir}/plink.${i}\_${j}.ge${no}.txt
        
            rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
            rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc
            rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.log
            rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
        
                fi
        fi

done < "$filename"

#merge effect alleles and AF to output files

echo "merge effect alleles and AF to output files"
cd ${section_16_dir}
sed 's/  \+/\t/g' < plink.${i}\_${j}.ge${no}.txt | sed 's/^ //g'  >plink.${i}.${PBS_ARRAYID}.ge${no}.txt
perl ${home_directory}/resources/phase2/join_file.pl -i "${section_16_dir}/plink.${i}.${PBS_ARRAYID}.ge${no}.txt,TAB,2 data.frq.tmp,TAB,1" -o ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt -a 1
#CpG    CHR SNP BP  NMISS   BETA    SE  R2  T   P   CHR SNP EA  NEA EAF N
awk '{print $1,$2,$3,$4,$13,$14,$15,$6,$7,$10,$16}' <${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt >${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt.tmp
mv ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt.tmp ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt
gzip ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt
rm plink.${i}\_${j}.ge${no}.txt plink.${i}.${PBS_ARRAYID}.ge${no}.txt
date

#Rscript ${home_directory}/resources/genetics/plotplinkvsgcta.R \
#    ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt.gz \
#    ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt.gz \
#    $j \
#    ${section_16_dir}
    
echo "Successfully completed"

