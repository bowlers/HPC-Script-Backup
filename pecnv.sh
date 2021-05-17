#!/bin/bash
#SBATCH --job-name=pecnv
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=kill-shared
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --error=./logs/pecnv/pecnv-%A-%a.err
#SBATCH --output=./logs/pecnv/pecnv-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load tools/matlab

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
mkdir "${PID}"
NORM_COUNTS="${DIR}"/baf/"${PID}"_normal_reads.txt
NORM_BAF="${DIR}"/baf/"${PID}"_normal.bed
TUM_COUNTS="${DIR}"/baf/"${PID}"_tumor_reads.txt
TUM_BAF="${DIR}"/baf/"${PID}"_tumor.bed

printf '%s\t%s\n' "${NORM_COUNTS}" "${NORM_BAF}" >> "${DIR}"/baf/"${PID}"/namelist.txt
printf '%s\t%s\n' "${TUM_COUNTS}" "${TUM_BAF}" >> "${DIR}"/baf/"${PID}"/namelist.txt
NAMELIST="${DIR}"/baf/"${PID}"/namelist.txt

OUTPUT="${DIR}"/baf/"${PID}"
GC="${RES}"/hg38.gc.bw
MAP="${RES}"/hg38.map.wig

matlab -nodisplay -nosplash -r 'try; cd("/home/sbowler/biocore_lts/scott/PECNV"); "workflow($NAMELIST,$OUTPUT,$GC,$MAP)"; catch; end; quit'
