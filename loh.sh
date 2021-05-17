#!/bin/bash
#SBATCH --job-name=LOH
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/loh/loh-%A-%a.err
#SBATCH --output=./logs/loh/loh-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/R

#R library definitions
if [ -n $R_LIBS ]; then
        export R_LIBS=/home/sbowler/biocore_lts/scott/R:$R_LIBS
else
        export R_LIBS=/home/sbowler/biocore_lts/scott/R
fi

LOC=/home/sbowler/biocore_lts/scott/R/CNV_Radar
RES=/home/sbowler/biocore_lts/scott/resources

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "Must provide a study"
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

if [ ! -d "${DIR}"/loh ]; then
	mkdir "${DIR}"/loh
fi

FILE="${DIR}"
FILE+=/design.txt

if [ ! -f "${DIR}"/loh/bam2roi.log ]; then
	if [ ! -f "${DIR}"/loh/normals.txt ]; then
		cut -f2 "${FILE}" >> "${DIR}"/loh/normals.txt
	fi
	while read p; do
		Rscript "${LOC}"/bam2roi.r -b "${DIR}"/bam/"${p}".bam -d "${RES}"/hg38.bed >> "${DIR}"/loh/bam2roi.log 2>&1
	done < "${DIR}"/loh/normals.txt
fi

