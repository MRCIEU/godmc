#!/bin/bash

set -e
source config


#Generate plink.raw / plink.frq #freq files are required to determine effect allele
echo "Converting plink files to transposed raw format"
${plink} --bfile ${bfile} --recode A-transpose --out ${bfile} --freq

echo "Generating matrixeqtl format genetic data"

# Getting the allele references for the 012 coding
cut -f 2,5,6 ${bfile}.traw > ${allele_ref}

# Get the header in the right format - just extract IID from the current header
head -n 1 ${bfile}.traw | cut -f 7- | tr '\t' '\n' | tr '_' '\t' | cut -f 2 | tr '\n' '\t' | awk '{ printf ("snpid\t%s\n", $0) }' > traw.header

# Remove extraneous columns
sed 1d ${bfile}.traw | cut -f 2,7- > ${bfile}.traw2
mv ${bfile}.traw2 ${bfile}.traw

echo "Splitting data into ${bfile_chunksize} SNP subsets"
split -d -a 10 -l ${bfile_chunksize} ${bfile}.traw ${tabfile}.tab.

i=1
for file in ${tabfile}.tab.*
do
	echo ${file}
	cat traw.header ${file} > ${tabfile}.tab.${i}
	rm ${file}
	i=$(($i+1))
done

rm traw.header


echo "Done!"
