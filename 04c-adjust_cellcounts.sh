set -e
source config

batch_number=${1}

if [ "${cellcounts_required}" = "no" ]
then
	echo "You have specified that cell count adjustment is not required. Writing appropriate files..."
	cp -v ${methylation_rt}.RData ${methylation_rt_cc}.RData
	echo "Done."
	exit
fi


if [ -n ${1} ]
then
	re='^[0-9]+$'
	if ! [[ $batch_number =~ $re ]] ; then
		echo "error: Batch variable is not a number"
		exit 1
	fi
	i=${1}
else
	i="NA"
fi

# Adjust for cell counts and rank transform
R --no-save --args ${betas} ${cellcounts} ${methylation_rt_cc} ${nthreads} ${meth_chunks} ${i} < resources/cellcounts/adjust_cellcounts.R
