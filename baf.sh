#!/bin/bash
#SBATCH --job-name=baf
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=kill-shared
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --error=./logs/baf/baf-%A-%a.err
#SBATCH --output=./logs/baf/baf-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/SAMtools
module load lang/Perl

RES=/home/sbowler/biocore_lts/scott/resources
if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "You must provide a study number."
	exit 1
else
        STUDY="${1}"
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"

if [ ! -d "${DIR}" ]; then
        echo "Invalid study.  Not located in /home/sbowler/biocore_lts/STUDY"
        exit 1
fi

BAND="/home/sbowler/biocore_lts/scott/resources/hg38.cytoband.txt"
FILE="${DIR}"
FILE+="/"
FILE+=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

if [ ! -d "${DIR}"/baf ]; then
	mkdir "${DIR}"/baf
fi

cd "${DIR}"/baf
if [ ! -f "${DIR}"/baf/"${PID}"_tumor.bed ]; then
	samtools mpileup -q 15 -f "${RES}"/GRCh38.fa -l "${RES}"/hg38.bed "${DIR}"/bam/"${TUMOR}".bam | ../../scripts/mpileup2baf.pl --min-reads 20 > "${PID}"_tumor.bed
	samtools view -c "${DIR}"/bam/"${TUMOR}".bam > "${PID}"_tumor_reads.txt
fi

if [ ! -f "${DIR}"/baf/"${PID}"_normal.ed ]; then
	samtools mpileup -q 15 -f "${RES}"/GRCh38.fa -l "${RES}"/hg38.bed "${DIR}"/bam/"${NORMAL}".bam | ../../scripts/mpileup2baf.pl --min-reads 20 > "${PID}"_normal.bed
	samtools view -c "${DIR}"/bam/"${NORMAL}".bam > "${PID}"_normal_reads.txt
fi

