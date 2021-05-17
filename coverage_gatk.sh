#!/bin/bash
#SBATCH --job-name=coverage_g
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --error=./logs/coverg-%A-%a.err
#SBATCH --output=./logs/coverg-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/env

if [ -z "${1}" ]; then
	echo You must provide a study ID. >&2
	exit 1
fi
if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi
if [ -z "${2}" ]; then
	TYPE=DNA
else
	TYPE=RNA
fi

DIR=/home/sbowler/biocore_lts/scott/
DIR+="${1}"
if [ ! -d "${DIR}" ]; then
	echo Could not locate directory. >&2
	exit 1
fi

RES=/home/sbowler/biocore_lts/scott/resources

if [ ! -d "${DIR}"/GATK ]; then
	mkdir "${DIR}"/GATK
fi
if [ ! -d "${DIR}"/GATK/coverage ]; then
	mkdir "${DIR}"/GATK/coverage
fi

cd "${DIR}"/bam

#  Get root filename
for i in `ls *.bam | sed -n 's/\.bam$//p'`; do
        FNAME="${i}"

	gatk DepthOfCoverage -R "${RES}"/hg38.fa -O ../GATK/coverage/"${FNAME}" -I "${FNAME}".bam -L "${RES}"/hg38.interval_list
done
