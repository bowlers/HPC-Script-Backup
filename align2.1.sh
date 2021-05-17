#!/bin/bash
#SBATCH --job-name=align2.1
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=56G
#SBATCH --error=./logs/align/align2.1-%A-%a.err
#SBATCH --output=./logs/align/align2.1-%A-%a.out


module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/Bowtie2
module load bio/SAMtools
module load lang/Python

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

if [ ! -d "${DIR}"/bam/metrics ]; then
	mkdir "${DIR}"/bam/metrics
fi

echo "${STUDY}"
echo "${FNAME}"

if [ ! -f "${DIR}"/bam/"${FNAME}".bam ]; then
	gatk MarkDuplicates -I "${DIR}"/bam/"${FNAME}".raw.bam -O "${DIR}"/bam/"${FNAME}".marked.bam -M "${DIR}"/bam/metrics/"${FNAME}".metrics.txt
	gatk BaseRecalibrator -I "${DIR}"/bam/"${FNAME}".marked.bam -R "${REF}" --known-sites $HOME/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz -O "${DIR}"/bam/"${FNAME}"_recal_data.table
	gatk ApplyBQSR -R "${REF}" -I "${DIR}"/bam/"${FNAME}".marked.bam --bqsr-recal-file "${DIR}"/bam/"${FNAME}"_recal_data.table -O "${DIR}"/bam/"${FNAME}".bam
	rm "${DIR}"/bam/"${FNAME}".raw.bam
	rm "${DIR}"/bam/"${FNAME}".marked.bam
	rm "${DIR}"/bam/"${FNAME}"_recal_data.table
	samtools index "${DIR}"/bam/"${FNAME}".bam
fi
