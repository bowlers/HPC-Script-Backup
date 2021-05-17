#!/bin/bash
#SBATCH --job-name=msi
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=shared
#SBATCH --time=00-01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --error=./logs/msi/msi-%A-%a.err
#SBATCH --output=./logs/msi/msi-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/SAMtools

if [ -z "${1}" ]; then
	echo You must provide a study ID.
	exit 1
else
	STUDY="${1}"
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

DIR=/home/sbowler/biocore_lts/scott/
DIR+="${1}"

if [ ! -d "${DIR}" ]; then
	echo Could not locate directory.
	exit 1
else
	cd "${DIR}"
fi

if [ ! -f ./design.txt ]; then
	echo Could not locate design.txt.
	exit 1
fi
DESIGN="${DIR}"/design.txt

AROW=$(sed -n "$SLURM_ARRAY_TASK_ID p" "${DESIGN}")
TUMOR=$( echo "${AROW}" | awk '{print $1}' )
NORMAL=$( echo "${AROW}" | awk '{print $2}' )
PID=$( echo "${AROW}" | awk '{print $3}' )

if [ ! -d ./msisensor ]; then
	mkdir "${DIR}"/msisensor
fi

if [ ! -d ./msisensor/"${PID}" ]; then
	mkdir "${DIR}"/msisensor/"${PID}"
fi
RES=/home/sbowler/biocore_lts/scott/resources

if [ ! -f "${RES}"/microsatellites.list ]; then
	msisensor scan -d "${RES}"/GRCh38.fa -o "${RES}"/microsatellites.list
fi

if [ ! -f ./bam/"${TUMOR}".bam.bai ]; then
	samtools index -@ 3 "${BAM}"/"${TUMOR}".bam
fi

if [ ! -f ./bam/"${NORMAL}".bam.bai ]; then
	samtools index -@ 3 "${BAM}"/"${NORMAL}".bam
fi

echo "${STUDY}"
echo "${PID}"

msisensor msi -b 2 -d "${RES}"/microsatellites.list -n ./bam/"${NORMAL}".bam -t ./bam/"${TUMOR}".bam -e "${RES}"/hg38.bed -o ./msisensor/"${PID}"/output
