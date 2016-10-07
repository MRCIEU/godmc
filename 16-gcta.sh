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

mydir="${methylation_processed_dir}"
cd $mydir

myout="${section_16_dir}"
echo $myout

i="1e-05"

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

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

nprobes=`awk -F' ' '{print NF; exit}'< ${methylation_processed_dir}/gcta.methylation.subset.${i}\_${j}.ge${no}.txt`
nprobes=$((nprobes - 2))

echo "$nprobes/$noprobes probes found in methylation data"
#echo "CpG" "CHR" "SNP" "BP" "NMISS" "BETA" "SE" "R2" "T" "P" | perl -pe 's/ /\t/g' > ${section_16_dir}/plink.${i}\_${j}.ge${no}.txt

Counter=0

while read -r line
do
    probe="$line"
    echo $probe

        probes=`grep -w $probe ${methylation_processed_dir}/gcta.methylation.subset.${i}\_${j}.ge${no}.txt | wc -l`

        if [ "$probes" -gt "0" ]   
        then
        
    zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $2}' > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
    
    #check whether fgrep is grepping one probe only
    #zcat ${home_directory}/resources/phase2/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt.gz | fgrep $probe | awk '{print $3}' |sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.probe
    
    awk -F":" '{print $1}' < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps | sort -u > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

    Counter=`expr $Counter + 1`
    echo "$Counter/$nprobes probes"



            while read -r line
            do
                chr="$line"
                echo $chr
      
                grep -w $chr ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps > ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
                nosnps=`grep -f ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr ${bfile}.bim | wc -l`
                echo "no snps = $nosnps"

                
                #check whether right probes are selected from methylation matrix
                probecheck=`expr $Counter + 2`
                probecheck2=`head -n1 ${methylation_processed_dir}/gcta.methylation.subset.${i}\_${j}.ge${no}.txt | cut -f $probecheck -d" "` 

                echo "$probe is $probecheck2"
                
                if [ "$nosnps" -gt "0" ]
                then
                chrno=`echo $chr |sed 's/chr//g'`     

                ${plink} \
                    --bfile ${bfile} \
                    --extract ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --make-bed \
                    --out ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
        
                if ! [[ "$probe" -eq "$probecheck2" ]] ; then
                echo "error: $probe doesn't match to methylation $probecheck2"
                exit 1
                fi

                
                if [ "$chrno" -lt "23" ]
                then

                ${gcta} \
                    --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --mlma \
                    --reml-no-constrain \
                    --grm ${grmfile_all}_minus_chr${chrno} \
                    --pheno ${methylation_processed_dir}/gcta.methylation.subset.${i}\_${j}.ge${no}.txt \
                    --mpheno $Counter \
                    --qcovar ${covariates_combined}.gcta.numeric \
                    --covar ${covariates_combined}.gcta.factor \
                    --out ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno} \
                    --thread-num ${nthreads}
                


                #    if [ ! "$?" -eq "0" ]
                 #   then 
                  #  ${plink} \
                   #     --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    #    --assoc \
                     #   --pheno ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt \
                      #  --mpheno $Counter \
                       # --out ${section_16_dir}/plink.${i}\_${probe}.ge${no}.chr${chrno} \
                        #--threads ${nthreads}
                    
                    #sed 's/^/'$probe'/' <${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
                    #tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${section_16_dir}/plink.${i}\_${j}.ge${no}.txt
        
                    #rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
                    #rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc
                    #rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.log



                    #fi
                fi
                
                if [ "$chrno" -eq "23" ]
                then
                ${gcta} \
                    --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                    --mlma \
                    --reml-no-constrain \
                    --grm ${grmfile_all} \
                    --pheno ${methylation_processed_dir}/gcta.methylation.subset.${i}\_${j}.ge${no}.txt \
                    --mpheno $Counter \
                    --qcovar ${covariates_combined}.gcta.numeric \
                    --covar ${covariates_combined}.gcta.factor \
                    --out ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno} \
                    --thread-num ${nthreads}
                
                    #if [ ! "$?" -eq "0" ]
                    #then 
                    #${plink} \
                     #   --bfile ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr \
                     #   --assoc \
                     #   --pheno ${methylation_processed_dir}/plink.methylation.subset.${i}\_${j}.ge${no}.txt \
                     #   --mpheno $Counter \
                     #   --out ${section_16_dir}/plink.${i}\_${probe}.ge${no}.chr${chrno} \
                     #   --threads ${nthreads}
                    
                    #sed 's/^/'$probe'/' <${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc | perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
                    #tail -n +2 ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp >>${section_16_dir}/plink.${i}\_${j}.ge${no}.txt
        
                    #rm ${methylation_processed_dir}/plink.${i}\_${probe}.ge${no}.qassoc.tmp
                    #rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.qassoc
                    #rm ${section_16_dir}/plink.${i}\_${probe}.ge${no}.log


                    #fi
                
                fi        


                echo "format 1"
                cat ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma | sed 's/^/'$probe'\t/'| perl -pe 's/  \+/ /g'  >${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
                echo "format 2"
                tail -n +2 ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp >>${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt
                echo "format 3" 
                rm ${methylation_processed_dir}/gcta.${i}\_${probe}.ge${no}.mlma.tmp
                rm ${section_16_dir}/gcta.${i}\_${probe}.ge${no}.chr${chrno}.mlma
                rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps.$chr
                #rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
                #rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.probe
                #rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs
                rm ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr.bed
                rm ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr.bim
                rm ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr.fam
                rm ${bfile}.${i}\_${probe}.ge${no}.allcohorts.snps.$chr.log
                
                fi
            
            done < ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs

        rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.snps
        rm ${genetic_processed_dir}/cis_trans.${i}\_${probe}.ge${no}.allcohorts.chrs
 
            
        fi
        
done < "$filename"

Nsubj=`grep observations ${section_16_dir}/logs_c/log.txt${batch_number} |tail -n1 |cut -d " " -f 1` 

awk '{if (NR==1) print $0,"N"; else print $0, '$Nsubj'}' OFS='\t' ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt > ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt.tmp
mv ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt.tmp ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt

nprobes2=`tail -n +2 ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt |awk '{print $1}'|sort -u |wc -l`
zerosnps=`grep "no snps = 0"  ${section_16_dir}/logs_c/log.txt${batch_number} |wc -l`
nprobes2=$(($nprobes2 + $zerosnps))

gzip ${section_16_dir}/gcta.${i}\_${j}.ge${no}.txt


if [ $nprobes2 -eq $nprobes ]
    then
    echo "Successfully completed ${batch_number}"
else
    echo "Completed $nprobes2 / $nprobes"

    exit 1
fi

#merge effect alleles and AF to output files

#echo "merge effect alleles and AF to plink output files"
#cd ${section_16_dir}
#sed 's/  \+/\t/g' < plink.${i}\_${j}.ge${no}.txt | sed 's/^ //g'  >plink.${i}\_${j}.ge${no}.txt.tmp
#mv plink.${i}\_${j}.ge${no}.txt.tmp plink.${i}\_${j}.ge${no}.txt
#cp plink.${i}\_${j}.ge${no}.txt plink.${i}.${PBS_ARRAYID}.ge${no}.txt
#perl ${home_directory}/resources/phase2/join_file.pl -i "${section_16_dir}/plink.${i}.${PBS_ARRAYID}.ge${no}.txt,TAB,2 data.frq.tmp,TAB,1" -o ${section_16_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt -a 1
#gzip ${section_17_dir}/plink.${i}\_${j}.ge${no}.gwama.formatted.txt
#rm plink.${i}\_${j}.ge${no}.txt plink.${i}.${PBS_ARRAYID}.ge${no}.txt


