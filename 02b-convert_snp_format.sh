#!/bin/bash

set -e
source config
exec &> >(tee ${section_02b_logfile})

function make_tab_format {

	transposed_file=$1
	allele_references=$2
	chunksize=$3
	outfile=$4

	echo "Generating matrixeqtl format genetic data"

	# Getting the allele references for the 012 coding
	cut -f 2,5,6 ${transposed_file}.traw > ${allele_references}

	# Get the header in the right format - just extract IID from the current header
	head -n 1 ${transposed_file}.traw | cut -f 7- | tr '\t' '\n' | tr '_' '\t' | cut -f 2 | tr '\n' '\t' | awk '{ printf ("snpid\t%s\n", $0) }' > traw.header

	# Remove extraneous columns
	sed 1d ${transposed_file}.traw | cut -f 2,7- > ${transposed_file}.traw2
	mv ${transposed_file}.traw2 ${transposed_file}.traw

	rm -f ${outfile}.tab.*
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
	rm ${transposed_file}.traw
	echo "Done!"

}



#Generate plink.raw / plink.frq #freq files are required to determine effect allele
echo "Converting plink files to transposed raw format"
${plink} --bfile ${bfile} --recode A-transpose --out ${bfile} --freq

# How big will each chunk be
nrow=`wc -l ${bfile}.bim | awk '{ print $1 }'`
chunksize=$(($nrow / $genetic_chunks))
remainder=$(($nrow % $genetic_chunks))

if [ ! "${remainder}" == "0" ]
then
	chunksize=$(($chunksize + 1))
fi
echo "Splitting genetic data into ${genetic_chunks} chunks"
echo "Each chunk will contain ${chunksize} SNPs"

make_tab_format ${bfile} ${allele_ref} ${chunksize} ${tabfile}

# Convert CNV data
Rscript resources/genetics/cnv_tabfile.R ${cnvs} ${tabcnv} ${genetic_chunks}

echo "Successfully converted genetic data"