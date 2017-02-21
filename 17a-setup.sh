#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_17a_logfile})
print_version


# Extract the SNP-CpG lists

echo "Downloading list of putative associations"

sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
get list_17.tgz
get list_17.tgz.md5sum
EOF

echo "Checking download integrity"

md5sum -c list_17.tgz.md5sum

echo "Extracting"

mkdir -p ${phase2_list_17}
tar xzf lists_17.tgz -C ${phase2_list_17}
rm lists_17.tgz*


# Create cpg files

head -n 1 ${methylation_adjusted_pcs}.txt > ${cpgfile17}_head.txt

for i in {1..100}
do
	c1=`cat ${phase2_list_17}/cpglist_${i}.txt | wc -l`
	echo "${i} of 100 - Extracting ${c1} probes"
	fgrep -wf ${phase2_list_17}/cpglist_${i}.txt ${methylation_adjusted_pcs}.txt > ${cpgfile17}_temp${i}.txt
	cat ${cpgfile17}_head.txt ${cpgfile17}_temp${i}.txt > ${cpgfile17}_list17_${i}.txt

	c2=`cat ${cpgfile17}_temp${i}.txt | wc -l`
	echo "Found ${c2} probes."
	rm ${cpgfile17}_temp${i}.txt
done

rm ${cpgfile17}_head.txt
echo "Done!"


echo "Creating genetic files"

nbatch=`ls -l ${tabfile}.tab.* | wc -l`
rm -f ${snpfile17}_temp
touch ${snpfile17}_temp
for i in $(seq 1 $nbatch)
do
	echo "Extracting SNPs from batch ${i} of ${nbatch}"
	fgrep -wf ${phase2_list_17}/snplist.txt ${tabfile}.tab.${i} >> ${snpfile17}_temp
done
head -n 1 ${tabfile}.tab.1 > ${snpfile17}_head
cat ${snpfile17}_head ${snpfile17}_temp > ${snpfile17}

c1=`cat ${phase2_list_17}/snplist.txt | wc -l`
c2=`cat ${snpfile17}_temp | wc -l`
echo "Found ${c2} out of ${c1} SNPs"
rm ${snpfile17}_temp ${snpfile17}_head

echo "Successfully completed script 17a"

