#!/bin/bash

set -e
source ./config

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

for i in {1..17}; do sections[$(($i-1))]=`printf "%02d" ${i}`; done
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
	temp=`which sshpassadsf | wc -l`

	if [ ! "${temp}" = "0" ]
	then
		read -s -p "Enter SFTP password: " mypassword
		export SSHPASS=${mypassword}
		echo "Testing connection"
		sshpass -e sftp -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path}/${sftp_username} << !
bye
!
		echo "Connection established"
	else
		echo "sshpass is not installed."
		echo "The results will now be archived, once that is done they will be uploaded to the server"
	fi

	echo ""
	echo "Tarring results and log files"
	tar czf results/${study_name}_${1}.tgz config resources/parameters results/${1}
	echo "Successfully created results archives"
	echo "Generating md5 checksum"
	md5sum results/${study_name}_${1}.tgz > results/${study_name}_${1}.md5sum
	
	echo ""

if [ ! "${temp}" = "0" ]
then

export SSHPASS=${mypassword}
sshpass -e sftp -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path}/${sftp_username} << !
   dir
   put results/${study_name}_${1}.md5sum
   put results/${study_name}_$1.tgz
   bye
!

else

read -s -p "Ready to upload? Press enter to continue: " anykey

sftp ${sftp_username}@${sftp_address}:${sftp_path}/${sftp_username} <<EOF
dir
put results/${study_name}_${1}.md5sum
put results/${study_name}_$1.tgz
EOF

fi


fi

