#!/bin/bash
#SBATCH --job-name=fastqc
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=shared
#SBATCH --time=00-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --error=./logs/fastqc-%A-%a.err
#SBATCH --output=./logs/fastqc-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/env
module load bio/FastQC

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
	echo Could not locate study. >&2
	exit 1
fi

if [ "${TYPE}" = "DNA" ]; then
	cd "${DIR}"/fastq
	if [ ! -d "${DIR}"/qc ]; then
		mkdir "${DIR}"/qc
	fi
else
	cd "${DIR}"/RNA/fastq
	if [ ! -d "${DIR}"/RNA/qc ]; then
		mkdir "${DIR}"/RNA/qc
	fi
fi

AROW=1
for i in `ls *.fastq | sed -n 's/\.fastq$//p'`;
        do
                if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
                        i=${i#"$PREFIX"}
                        FNAME="${i}"
                        break
                fi
                ((AROW++))
        done

fastqc -o ../qc "${FNAME}".fastq
