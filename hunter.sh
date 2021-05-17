#!/bin/bash
#SBATCH --job-name=hunter
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/hunter/hunter-%A-%a.err
#SBATCH --output=./logs/hunter/hunter-%A-%a.out

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env
module load bio/SAMtools
module load lang/R


#R library definitions
if [ -n $R_LIBS ]; then
        export R_LIBS=/home/sbowler/biocore_lts/scott/Rlib:$R_LIBS
else
        export R_LIBS=/home/sbowler/biocore_lts/scott/Rlib
fi

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

if [ ! -d "${DIR}"hunter ]; then
	mkdir "${DIR}"/hunter
fi

cd "${DIR}"

BAND="$HOME/biocore_lts/scott/resources/hg38.cytoband.txt"
FILE=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

cd ./bam

if [ ! -f ./bam/"${TUMOR}".bam.bai ]; then
	samtools index -b "${TUMOR}".bam
fi
if [ ! -f ./bam/"${NORMAL}".bam.bai ]; then
	samtools index -b "${NORMAL}".bam
fi
cd ..

echo "${STUDY}"
echo "${PID}"

if [ ! -f ./hunter/"${PID}"/"${PID}"_summary.tsv ]; then
	telomerehunter -b "${BAND}" -ibt ./bam/"${TUMOR}".bam -ibc ./bam/"${NORMAL}".bam -o ./hunter -p "${PID}" -p1
fi
