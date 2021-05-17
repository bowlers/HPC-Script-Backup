#!/bin/bash
#SBATCH --job-name=TCR
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/tcr/tcr-%A-%a.err
#SBATCH --output=./logs/tcr/tcr-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env
module load bio/SAMtools
module load lang/Java
export PATH="$HOME/biocore_lts/scott/mixcr:$PATH"

if [ -z "${1}" ]; then
        echo You must provide a study ID. >&2
        exit 1
fi
STUDY="${1}"

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

DIR=/home/sbowler/biocore_lts/scott/
DIR+="${1}"
if [ ! -d "${DIR}" ]; then
        echo Could not locate directory. >&2
        exit 1
fi

RES=/home/sbowler/biocore_lts/scott/resources
REF="${RES}"/GRCh38.fa
cd "${DIR}"

#  Get root filename
PREFIX="./fastq/"
AROW=1

for i in `ls ./fastq/*_1.fastq | sed -n 's/\_1.fastq$//p'`; do
        if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
                                i=${i#"$PREFIX"}
	                        FNAME="${i}"
                                break
        fi
        ((AROW++))
done

if [ ! -d "${DIR}"/tcr ]; then
        mkdir "${DIR}"/tcr
fi

echo "${STUDY}"
echo "${FNAME}"

mixcr analyze amplicon --species hs --starting-material dna --5-end v-primers --3-end j-primers --adapters no-adapters ./fastq/"${FNAME}"_1.fastq ./fastq/"${FNAME}"_2.fastq ./tcr/"${FNAME}"

