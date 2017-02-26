#!/bin/bash

set -e
source ./config
exec &> >(tee ${section_17a_logfile})
print_version


# Extract the SNP-CpG lists

echo "Downloading list of putative associations"

sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
get lists_17.tgz
get lists_17.tgz.md5sum
EOF

echo "Checking download integrity"

md5sum -c lists_17.tgz.md5sum

echo "Extracting"

rm -rf ${phase2_list_17}
mkdir -p ${phase2_list_17}
tar xzf lists_17.tgz -C ${phase2_list_17}
rm lists_17.tgz*




echo "Creating genetic files"

n_genetic_batch=`ls -l ${tabfile}.tab.* | wc -l`
echo "Genetic data is split into ${n_genetic_batch} chunks"

if [ ! "${n_genetic_batch}" = "${genetic_chunks}" ]
then
	echo "Problem: Genetic data has been split into ${n_genetic_batch}, but the number of batches specified in the config file is ${genetic_chunks}"
	echo "Please either change the 'genetic_chunks' variable in the config file to ${n_genetic_batch} or re-run script 02b"
	exit
fi


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


echo "Creating SNP info file"

${plink} --bfile ${bfile} \
	--freq gz \
	--extract ${phase2_list_17}/snplist.txt \
	--out ${phase2_list_17}/snpinfo


echo "Successfully completed script 17a"

