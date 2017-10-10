#!/bin/bash

set -e
source ./config


        if [ -f "processed_data/genetic_data/data.easyqc.flipped.SNPs.txt" ]; then
                echo "easy qc flipped snps file present"
        else
            	echo "Problem: easy qc flipped snps file is absent"
                exit 1
        fi

        nsnps=`wc -l processed_data/genetic_data/data.easyqc.flipped.SNPs.txt |awk '{print $1}'`
	echo $nsnps present in processed_data/genetic_data/data.easyqc.flipped.SNPs.txt
                if [ "$nsnps" -gt "0" ]; then
                        echo "easy qc flipped snps file is not empty"
                else
                    	echo "easy qc flipped snps file is empty"
                        exit 1
                fi


	echo ""
	temp=`which sshpass 2>/dev/null | wc -l`

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

	echo ""

if [ ! "${temp}" = "0" ]
then

export SSHPASS=${mypassword}
sshpass -e sftp -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path}/${sftp_username} << !
   dir
   put processed_data/genetic_data/data.easyqc.flipped.SNPs.txt   
   bye
!

else

read -s -p "Ready to upload? Press enter to continue: " anykey

sftp ${sftp_username}@${sftp_address}:${sftp_path}/${sftp_username} <<EOF
dir
put processed_data/genetic_data/data.easyqc.flipped.SNPs.txt
EOF

fi


