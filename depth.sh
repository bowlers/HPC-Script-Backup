#!/bin/bash
#SBATCH --job-name=depth
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/depth/depth-%A-%a.err
#SBATCH --output=./logs/depth/depth-%A-%a.out

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env
module load lang/Python

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "Provide a study"
        exit 1
else
        STUDY="${1}"
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"

if [ ! -d "${DIR}" ]; then
	echo "Coult not locate study directory"
	exit 1
else
	cd "${DIR}"
fi

if [ ! -d ./depth ]; then
	mkdir depth
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" design.txt)
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

INTERVAL=$HOME/biocore_lts/scott/resources/hg38_WES.interval_list
REF=$HOME/biocore_lts/scott/resources/GRCh38.fa

gatk DepthOfCoverage --input ./bam/"${TUMOR}".bam --intervals "${INTERVAL}" --output ./depth/"${TUMOR}" --reference "${REF}"
gatk DepthofCoverage --input ./bam/"${NORMAL}".bam -- intervals "${INTERVAL}" --output ./depth/"${NORMAL}" --reference "${REF}"
