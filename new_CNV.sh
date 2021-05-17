#!/bin/bash
#SBATCH --job-name=CNV_ADTEx
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --error=./logs/ADTEx-%A-%a.err
#SBATCH --output=./logs/ADTEx-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/R/3.5.1-intel-2018.5.274-Python-2.7.15
module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/env
module load lang/Python/2.7.15-GCCcore-8.2.0

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
        STUDY=phs000452
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

FILE="${DIR}"
FILE+=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

cd $HOME/biocore_lts/scott/ADTEx.v.2.0
python ADTEx.py --normal "${DIR}"bam/"${NORMAL}".bam --tumor "${DIR}"bam/"${TUMOR}".bam --bed $HOME/biocore_lts/scott/resources/hg38.cytoband.bed --out "${DIR}"bedtools/"${PID}"
