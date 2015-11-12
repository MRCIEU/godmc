#!/bin/bash

set -e
source config

checkFirstArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid section identifier"
	echo $"Need to specify a value from 01 to $e"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}

checkSecondArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $2 is not a valid action"
	echo $"Specify either 'check' or 'upload'"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}

source resources/logs/check_logs.sh
source resources/logs/check_results.sh

for i in {1..13}; do sections[$(($i-1))]=`printf "%02d" ${i}`; done
checkFirstArg "$1" "${sections[@]}"

actions=("check" "upload")
checkSecondArg "$2" "${actions[@]}"

echo ""
echo "Checking log files for $1"
eval "check_logs_$1"

echo ""
echo "Checking results for $1"
eval "check_results_$1"

echo ""
echo "Section $1 has been successfully completed!"

if [ "$2" = "upload" ]
then
	echo ""
	echo "Tarring results and log files"
	tar czf results/uploads_$1.tgz config resources/parameters results/$1
	echo "Successfully created results archives"
	echo "Generating md5 checksum"
	md5sum results/uploads_$1.tgz > results/md5sum_$1.txt
	
	echo ""

sftp ${sftp_username}@${sftp_address}:/${sftp_username} <<EOF
dir
put results/md5sum_$1.txt
put results/uploads_$1.tgz
EOF

fi

