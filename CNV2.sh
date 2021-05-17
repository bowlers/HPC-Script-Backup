#!/bin/bash
#SBATCH --job-name=CNV_smooth
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --error=./logs/cnv/cnv2-%A-%a.err
#SBATCH --output=./logs/cnv/cnv2-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/R

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
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

if [ ! -d "${STUDY_DIR}"/loh ]; then
	echo "LoH has not been performed"
	exit 1
fi

FILE="${STUDY_DIR}"
FILE+=/design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

Rscript CNV.R "${STUDY}" "${PID}"
