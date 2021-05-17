#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --error=./logs/bam/cover-%A-%a.err
#SBATCH --output=./logs/bam/cover-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env

RES=$HOME/biocore_lts/scott/resources/hg38.bed

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

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ ! -d "${STUDY_DIR}"/cover ]; then
        mkdir "${STUDY_DIR}"/cover
fi

cd "${STUDY_DIR}"
FILE="${DIR}"
FILE+=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

echo "${STUDY}"
echo "${PID}"

if [ ! -f "${STUDY_DIR}"/cover/"${TUMOR}".coverage.tsv ]; then
	samtools depth -H -o "${STUDY_DIR}"/cover/"${TUMOR}".coverage.tsv "${STUDY_DIR}"/bam/"${TUMOR}".bam
fi
