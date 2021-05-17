#!/bin/bash
#SBATCH --job-name=egaUpdate
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
##SBATCH --partition=shared
#SBATCH --time=02-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/ega/egaUp-%A-%a.err
#SBATCH --output=./logs/ega/egaUp-%A-%a.out

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

if [ "${2}" == "rna" ]; then
	if [ ! -d ./RNA ]; then
		echo Could not locate RNA folder.
		exit 1
	fi
	cd RNA
fi

if [ ! -d ./bam ]; then
	echo Could not locate bam folder.
	exit 1
fi

echo $STUDY
pwd

for i in `ls *.*`; do

	if [ "$i" == "EGA_update_names.sh" ]; then
		continue
	fi

	EXT=${i: -8}

	if [ "$EXT" != ".bam.bai" ]; then
		EXT=${i: -4}
	fi

	NEW=$(echo $i | sed -r 's/_/\t/g')
	NEW2=$(echo $NEW | awk '{print $(NF-1)}')
	NEW2+=$EXT
	mv $i $NEW2
done
