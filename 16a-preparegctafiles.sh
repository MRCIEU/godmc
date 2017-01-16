#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_16a_logfile})
print_version

#pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")

# 74 subgroups
#cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

# 128 subgroups
cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-2]" "cg01[3-5]" "cg01[6-8]" "cg01[9]" "cg02[0-2]" "cg02[3-5]" "cg02[6-8]" "cg02[9]" "cg03[0-2]" "cg03[3-5]" "cg03[6-8]" "cg03[9]" "cg04[0-2]" "cg04[3-5]" "cg04[6-8]" "cg04[9]" "cg05[0-2]" "cg05[3-5]" "cg05[6-8]" "cg05[9]" "cg06[0-2]" "cg06[3-5]" "cg06[6-8]" "cg06[9]" "cg07[0-2]" "cg07[3-5]" "cg07[6-8]" "cg07[9]" "cg08[0-2]" "cg08[3-5]" "cg08[6-8]" "cg08[9]" "cg09[0-2]" "cg09[3-5]" "cg09[6-8]" "cg09[9]" "cg10[0-2]" "cg10[3-5]" "cg10[6-8]" "cg10[9]" "cg11[0-2]" "cg11[3-5]" "cg11[6-8]" "cg11[9]" "cg12[0-2]" "cg12[3-5]" "cg12[6-8]" "cg12[9]" "cg13[0-2]" "cg13[3-5]" "cg13[6-8]" "cg13[9]" "cg14[0-2]" "cg14[3-5]" "cg14[6-8]" "cg14[9]" "cg15[0-2]" "cg15[3-5]" "cg15[6-8]" "cg15[9]" "cg16[0-2]" "cg16[3-5]" "cg16[6-8]" "cg16[9]" "cg17[0-2]" "cg17[3-5]" "cg17[6-8]" "cg17[9]" "cg18[0-2]" "cg18[3-5]" "cg18[6-8]" "cg18[9]" "cg19[0-2]" "cg19[3-5]" "cg19[6-8]" "cg19[9]" "cg20[0-2]" "cg20[3-5]" "cg20[6-8]" "cg20[9]" "cg21[0-2]" "cg21[3-5]" "cg21[6-8]" "cg21[9]" "cg22[0-2]" "cg22[3-5]" "cg22[6-8]" "cg22[9]" "cg23[0-2]" "cg23[3-5]" "cg23[6-8]" "cg23[9]" "cg24[0-2]" "cg24[3-5]" "cg24[6-8]" "cg24[9]" "cg25[0-2]" "cg25[3-5]" "cg25[6-8]" "cg25[9]" "cg26[0-2]" "cg26[3-5]" "cg26[6-8]" "cg26[9]" "cg27[0-2]" "cg27[3-5]" "cg27[6-8]" "cg27[9]" "_ch")

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



#for i in `seq 1 22`; do
#
#echo "Creating kinship matrices minus chromosome ${i}"
#zcat ${hm3_snps} | grep -v -w chr${i} |sort -u > $genetic_processed_dir/grm_minus_chr${i}.snps
#
#${plink} \
#	--bfile ${bfile} \
#	--extract $genetic_processed_dir/grm_minus_chr${i}.snps \
#	--maf ${grm_maf_cutoff} \
#	--make-grm-bin \
#	--out ${grmfile_all}_minus_chr${i} \
#	--threads ${nthreads} \
#	--autosome
#
#rm $genetic_processed_dir/grm_minus_chr${i}.snps
#
#done
#
#echo "Creating files for gcta"
#
#	Rscript ./resources/methylation/extractprobesets3.R \
#	${betas} \
#	${covariates_combined}.txt \
#	${nongenetic_meth_pcs}.txt \
#	${#cpgs[@]} \
#	./resources/phase2 \
#	${methylation_processed_dir} \
#	${bfile}.fam \
#	${covariates_combined}.gcta \
#	${nthreads}

echo "Creating files for plink"
Rscript ./resources/methylation/convertbetamatrix.R \
	${methylation_adjusted_pcs}.RData \
	${#cpgs[@]} \
	${methylation_processed_dir} \
	./resources/phase2 \
	${bfile}.fam

echo "Calculating MAF"
${plink} \
    --bfile ${bfile} \
    --freq gz \
    --out ${section_16_dir}/data

zcat ${section_16_dir}/data.frq.gz | sed -e 's/[[:space:]]\+/ /g' |perl -pe 's/^ //g'|perl -pe 's/ /\t/g'|awk -v OFS='\t' '{ if(NR>1) print $1,$2,$3,$4,$5,$6/2; else print $0;}'|perl -pe 's/A1/EA/g' |perl -pe 's/A2/NEA/g' |perl -pe 's/MAF/EAF/g'|perl -pe 's/NCHROBS/N/g' |perl -pe 's/ /\t/g'>${section_16_dir}/data.frq.tmp

echo "Successfully completed script 16a"

