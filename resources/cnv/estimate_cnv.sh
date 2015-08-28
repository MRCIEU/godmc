#!/bin/bash


# The number of jobs will be based on 
# how many batches you want to split the
# samples across. 
# 
# For example, for 500 samples, if I want
# each node to run 10 samples and using
# 2 cores per node then I would specify 
# the -t flag to be 1-50, the ppn value
# to be 2, and the $splitsize variable to 
# be 10.

 
#PBS -N cnv
#PBS -o cnv-output
#PBS -e cnv-error
#PBS -t 1-50
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
 
set -e
 
echo "Running on ${HOSTNAME}"
 
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi
 
i=${PBS_ARRAYID}


#####################
# For users to edit #
#####################

# Path to godmc directory
godmc_dir=""

# Path to idat files
idat_dir="/path/to/idat_files"

# Number of cores to use per batch
ncores="2"

# Number of samples to calculate per batch
splitsize="10"

# Output file rootname
output_name="cnv_filename"


#######################################
# No editing required from this point #
#######################################

output="${godmc_dir}/raw_data/cnv/${output_name}_${i}.RData"
script="${godmc_dir}/resources/cnv/estimate_cnv.R"


R --no-save --args ${idat_dir} ${output} ${ncores} ${i} ${splitsize} < ${script}

