#!/bin/bash
#SBATCH --job-name=move
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/move/move-%A-%a.err
#SBATCH --output=./logs/move/move-%A-%a.out

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "You must provide a study ID."
        exit 1
else
        STUDY="${1}"
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"
DIR+="/"

if [ ! -d "${DIR}" ]; then
        echo "Invalid study.  Not located in /home/sbowler/biocore_lts/STUDY"
        exit 1
fi

cd "${DIR}"

if [ ! -d ./fastq/ ]; then
	mkdir fastq
fi

if [ ! -f files.txt ]; then
	echo "Could not locate files list."
	exit 1
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" files.txt)

if [ -d ./"${LINE}" ]; then
	cd "${LINE}"

	if ls ./*.gz &> /dev/null; then
		gzip -d *.gz
		mv *.fastq ../fastq/
	fi
	if ls ./*.md5 &> /dev/null; then
		if [ ! -d ../md5 ]; then
			mkdir ../md5
		fi
		mv *.md5 ../md5
	fi
	cd ..
	rm -r "${LINE}"
fi
