#!/bin/bash

set -e
source config

function make_tab_format {

	plink_bin=$1
	transposed_file=$2
	allele_references=$3
	chunksize=$4
	outfile=$5

	echo "Generating matrixeqtl format genetic data"

	# Getting the allele references for the 012 coding
	cut -f 2,5,6 ${transposed_file}.traw > ${allele_references}

	# Get the header in the right format - just extract IID from the current header
	head -n 1 ${transposed_file}.traw | cut -f 7- | tr '\t' '\n' | tr '_' '\t' | cut -f 2 | tr '\n' '\t' | awk '{ printf ("snpid\t%s\n", $0) }' > traw.header

	# Remove extraneous columns
	sed 1d ${transposed_file}.traw | cut -f 2,7- > ${transposed_file}.traw2
	mv ${transposed_file}.traw2 ${transposed_file}.traw

	echo "Splitting data into ${chunksize} SNP subsets"
	rm ${outfile}.tab.*
	split -d -a 10 -l ${chunksize} ${transposed_file}.traw ${outfile}.tab.

	i=1
	for file in ${outfile}.tab.*
	do
		echo ${file}
		cat traw.header ${file} > ${outfile}.tab.${i}
		rm ${file}
		i=$(($i+1))
	done

	rm traw.header
	echo "Done!"

}



#Generate plink.raw / plink.frq #freq files are required to determine effect allele
echo "Converting plink files to transposed raw format"
${plink} --bfile ${bfile} --recode A-transpose --out ${bfile} --freq
make_tab_format ${plink} ${bfile} ${allele_ref} ${bfile_chunksize} ${tabfile}


echo "Converting hapmap3 SNPs to transposed raw format"
gunzip -c ${hm3_snps} > temp.snplist
${plink} --bfile ${bfile} --recode A-transpose --out ${bfile_hm3} --extract temp.snplist
rm temp.snplist
make_tab_format ${plink} ${bfile_hm3} ${allele_ref_hm3} ${bfile_chunksize} ${tabfile_hm3}


