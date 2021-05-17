#!/bin/bash
#SBATCH --job-name=align2
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --error=./logs/align/align2-%A-%a.err
#SBATCH --output=./logs/align/align2-%A-%a.out


module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/Bowtie2
module load bio/SAMtools

if [ -z "${1}" ]; then
	echo You must provide a study ID. >&2
	exit 1
fi
STUDY="${1}"

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

DIR=/home/sbowler/biocore_lts/scott/
DIR+="${1}"
if [ ! -d "${DIR}" ]; then
	echo Could not locate directory. >&2
	exit 1
fi

RES=/home/sbowler/biocore_lts/scott/resources 
REF="${RES}"/GRCh38.fa 
cd "${DIR}"

#  Get root filename
PREFIX="./fastq/"
AROW=1

if [ -z "${2}" ]; then
	for i in `ls ./fastq/*_1.fastq | sed -n 's/\_1.fastq$//p'`;
        	do
                	if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
                        	i=${i#"$PREFIX"}
	                        FNAME="${i}"
        	                break
                	fi
	                ((AROW++))
        	done
else
	FNAME="${2}"
fi

if [ ! -d "${DIR}"/bam ]; then
	mkdir "${DIR}"/bam
fi

if [ ! -d "${DIR}"/bam/metrics ]; then
	mkdir "${DIR}"/bam/metrics
fi

echo "${STUDY}"
echo "${FNAME}"

if [ ! -f "${DIR}"/bam/"${FNAME}".bam ]; then
	cd "${RES}"
	bowtie2 -p 4 --rg-id 1 --rg LB:lib1 --rg PL:ILLUMINA --rg PU:unit1 --rg SM:"${FNAME}" -x GRCh38 -1 "${DIR}"/fastq/"${FNAME}"_1.fastq -2 "${DIR}"/fastq/"${FNAME}"_2.fastq -S "${DIR}"/bam/"${FNAME}".sam
	cd "${DIR}"/bam
	samtools view -bS "${FNAME}".sam > "${FNAME}"_unsorted.bam
	samtools sort -@ 4 "${FNAME}"_unsorted.bam -o "${FNAME}".raw.bam
	rm "${FNAME}"_unsorted.bam
	rm "${FNAME}".sam

	cd $HOME/biocore_lts/scott/scripts
	if [ -z "${2}" ]; then
		jid=$(sbatch --array="${SLURM_ARRAY_TASK_ID}" align2.1.sh "${STUDY}")
	else
		jid=$(sbatch align2.1.sh "${STUDY}" "${FNAME}")
	fi
fi
