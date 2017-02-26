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
exec &> >(tee ${section_17b_logfile}${batch_number})
print_version


i=${batch_number}

nbatch=(${phase2_list_17}/cpglist_*.txt)
nbatch=${#nbatch[@]}


# Create cpg files

head -n 1 ${methylation_adjusted_pcs}.txt > ${cpgfile17}_head_${i}.txt

c1=`cat ${phase2_list_17}/cpglist_${i}.txt | wc -l`
echo "Performing ${i} of ${nbatch}"
echo "Extracting ${c1} probes"
fgrep -wf ${phase2_list_17}/cpglist_${i}.txt ${methylation_adjusted_pcs}.txt > ${cpgfile17}_temp${i}.txt
cat ${cpgfile17}_head_${i}.txt ${cpgfile17}_temp${i}.txt > ${cpgfile17}_list17_${i}.txt

c2=`cat ${cpgfile17}_temp${i}.txt | wc -l`
echo "Found ${c2} probes."
rm ${cpgfile17}_temp${i}.txt ${cpgfile17}_head_${i}.txt


geno="${snpfile17}"
phen="${cpgfile17}_list17_${i}.txt"
cov="NULL"
threshold="1"
out="${section_17_dir}/results_${i}.txt.gz"


echo "Performing meQTL analysis"

Rscript resources/genetics/matrixeqtl_small.R \
	${geno} \
	${phen} \
	${cov} \
	${threshold} \
	${phase2_list_17}/snpinfo.frq.gz \
	${out}

rm ${phen}

echo "Successfully completed ${batch_number} of ${nbatch}"
