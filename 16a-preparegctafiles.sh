#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16a_logfile})
print_version

#pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

#pval cut-off
i="1e-05"

#cpgs=("cg0000[0-9]")

#found in number of cohorts: cis=1 trans=2
no="1.2"


	sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
mget cis_trans.1e-05_*.ge${no}.allcohorts.txt.gz
mget cis_trans.1e-05_*.ge${no}.allcohorts.probes
EOF

	mv cis_trans.1e-05_*.ge${no}.allcohorts.txt.gz ${home_directory}/resources/phase2
	mv cis_trans.1e-05_*.ge${no}.allcohorts.probes ${home_directory}/resources/phase2



for i in `seq 1 22`; do

echo "Creating kinship matrices minus chromosome ${i}"
zcat ${hm3_snps} | grep -v -w chr${i} |sort -u > $genetic_processed_dir/grm_minus_chr${i}.snps

${plink} \
	--bfile ${bfile} \
	--extract $genetic_processed_dir/grm_minus_chr${i}.snps \
	--maf ${grm_maf_cutoff} \
	--make-grm-bin \
	--out ${grmfile_all}_minus_chr${i} \
	--threads ${nthreads} \
	--autosome

#rm $genetic_processed_dir/grm_minus_chr${i}.snps

done

echo "Creating files for gcta"

Rscript ./resources/methylation/extractprobesets3.R \
	${betas} \
	${covariates_combined}.txt \
	${#cpgs[@]} \
	./resources/phase2 \
	${methylation_processed_dir} \
	${bfile}.fam \
	${covariates_combined}.gcta

#echo "Calculating MAF"
#${plink} \
#    --bfile ${bfile} \
#    --freq gz \
#    --out ${lmm_res_dir}/data

#zcat ${lmm_res_dir}/data.frq.gz | sed -e 's/[[:space:]]\+/ /g' |perl -pe 's/^ //g'|perl -pe 's/ /\t/g'|awk -v OFS='\t' '{ if(NR>1) print $1,$2,$3,$4,$5,$6/2; else print $0;}'|perl -pe 's/A1/EA/g' |perl -pe 's/A2/NEA/g' |perl -pe 's/MAF/EAF/g'|perl -pe 's/NCHROBS/N/g' |perl -pe 's/ /\t/g'>${lmm_res_dir}/data.frq.tmp


echo "Successfully completed script 16a"

