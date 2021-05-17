#!/bin/bash
#SBATCH --job-name=varscan
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --error=./logs/TMB/vTMB-%A-%a.err
#SBATCH --output=./logs/TMB/vTMB-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/SAMtools

RES=$HOME/biocore_lts/scott/resources/GRCh38.fa

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

STUDY="${1}"
if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
        echo "Could not locate study directory"
        exit 1
fi
STUDY_DIR=$HOME/biocore_lts/scott/"${STUDY}"

# Locate design file so that we can pull normal and tumor SRRs and PID
if [ ! -f "${STUDY_DIR}"/design.txt ]; then
        echo No design file within study. >&2
        exit 1
fi

if [ ! -d "${STUDY_DIR}"/varscan ]; then
	mkdir "${STUDY_DIR}"/varscan
fi


# Design exists, so lets pull out the tumor and normal SRRs and PID
DESIGN="${STUDY_DIR}/design.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${DESIGN}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

cd "${STUDY_DIR}"/varscan
samtools mpileup -q 1 -f "${RES}" "${STUDY_DIR}"/bam/"${NORMAL}".bam > "${NORMAL}".pileup
samtools mpileup -q 1 -f "${RES}" "${STUDY_DIR}"/bam/"${TUMOR}".bam > "${TUMOR}".pileup

varscan somatic "${STUDY_DIR}"/varscan/"${NORMAL}".pileup "${STUDY_DIR}"/varscan/"${TUMOR}".pileup "${STUDY_DIR}"/varscan/"${PID}" --somatic-p-value 0.05 --p-value 0.99 --output-vcf 1

rm "${NORMAL}".pileup
rm "${TUMOR}".pileup
