set -e
source config

batch_number=${1}

if [ "${unrelated}" = "yes" ]
then
	echo "You have specified that the data is not family data. Writing appropriate files..."
	cp -v ${methylation_rt_cc} ${methylation_rt_poly}
	echo "Done."
	exit
fi


if [ -n "${1}" ]
then
	re='^[0-9]+$'
	if ! [[ $batch_number =~ $re ]] ; then
		echo "error: Batch variable is not a number"
		exit 1
	fi
	i=${1}
	echo "Running batch ${i} of ${meth_chunks}"
else
	i="NA"
	echo "Running entire set on a single node using ${nthreads} threads."
fi



# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)

R --no-save --args ${methylation_rt_cc} ${grmfile_relateds} ${methylation_rt_cc_poly} ${nthreads} ${meth_chunks} ${i} < resources/relateds/methylation_relateds.R

