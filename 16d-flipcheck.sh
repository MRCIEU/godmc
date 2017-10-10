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
	cp processed_data/genetic_data/data.easyqc.flipped.SNPs.txt processed_data/genetic_data/${study_name}_data.easyqc.flipped.SNPs.txt

	echo ""

	        sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources <<EOF
put processed_data/genetic_data/${study_name}_data.easyqc.flipped.SNPs.txt
EOF




