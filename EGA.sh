#!/bin/bash
#SBATCH --job-name=egadownload
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=shared
#SBATCH --time=02-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/ega/ega-%A-%a.err
#SBATCH --output=./logs/ega/ega-%A-%a.out

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/python3.7

if [ -z "${1}" ]; then
	echo You must provide an EGAD number.
	exit 1
else
	STUDY="${1}"
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
	SLURM_ARRAY_TASK_ID=1
fi

cd $HOME/biocore_lts/scott

if [ ! -d ./"${STUDY}" ]; then
	echo Are you sture?  This study directory does not exist.
	exit 1
fi

cd "${STUDY}"

if [ ! -f ./files.txt ]; then
	echo No files.txt to list all known files.
	exit 1
fi

if [ -z "${2}" ] || [ "${2}" == "rna" ]; then
	FILE=files.txt
	FNAME=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")

	if [ ! -d ./"${FNAME}" ]; then
		echo "Could not locate directory, download file"
		pyega3 -cf ../credential.json fetch "${FNAME}"
	else
		COUNT=`ls -1 ./"${FNAME}"/*.md5 2>/dev/null | wc -l`
		if [ "${COUNT}" != 0 ]; then
			echo "MD5 found!"
			cd "${FNAME}"
			MD5=`ls *.md5`
			tMD5=${MD5%.gz.md5}
			cd ..
			if [ "${2}" == "rna" ]; then
	                        if [ ! -f ./RNA//fastq/"${tMD5}" ]; then
        	                        echo "FASTQ Not found, Download!"
                	                pyega3 -cf ../credential.json fetch "${FNAME}"
                        	else
                                	echo "FASTQ found, skipping."
	                        fi
			else
				if [ ! -f ./fastq/"${tMD5}" ]; then
					echo "FASTQ Not found, Download!"
					pyega3 -cf ../credential.json fetch "${FNAME}"
				else
					echo "FASTQ found, skipping."
				fi
			fi
		else
			echo "MD5 not found, download file!"
			pyega3 -cf ../credential.json fetch "${FNAME}"
		fi
	fi
#	pyega3 -cf ../credential.json fetch "${FNAME}"
else
	pyega3 -cf ../credential.json fetch "${2}"
fi
