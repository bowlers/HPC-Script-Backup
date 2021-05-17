#!/bin/bash
#SBATCH --job-name=mixcr
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/mixcr/mixcr-%A-%a.err
#SBATCH --output=./logs/mixcr/mixcr-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env

export PATH=$HOME/biocore_lts/scott/mixcr:$PATH

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "Syntax: squeue (optional: --array=x-y) mixcr.sh <study> <type>"
        echo "Array and type parameters are optional"
        echo "Valid type arguments: dna rna"
        exit 1
else
	STUDY="${1}"
fi


if [ -z "${2}" ]; then
	TYPE=dna
else
	TYPE="${2}"
	if [ $TYPE != "dna" ] && [ $TYPE != "rna" ]; then
		echo "Syntax: squeue (optional: --array=x-y) mixcr.sh <study> <type>"
		echo "Array and type parameters are optional"
		echo "Valid type arguments: dna rna"
		exit 1
	fi	
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"
DIR+="/"

if [ ! -d "${DIR}" ]; then
        echo "Invalid study.  Not located in /home/sbowler/biocore_lts/STUDY"
        exit 1
fi

FILE="${DIR}"
cd "${DIR}"

if [ $TYPE == "dna" ]; then
	FILE+=design.txt
	LOC="${DIR}"fastq
	if [ ! -d "${DIR}"/mixcr ]; then
		mkdir "${DIR}"/mixcr
	fi
	LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
	TUMOR=$(echo $LINE | awk '{print $1}')
	NORMAL=$(echo $LINE | awk '{print $2}')
	PID=$(echo $LINE | awk '{print $3}')
	OUTPUT=./mixcr/"${PID}"
else
	FILE+=RNA/design.txt
	LOC="${DIR}"RNA/fastq
	if [ ! -d "${DIR}"RNA/mixcr ]; then
		mkdir "${DIR}"RNA/mixcr
	fi
	LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")

	NORMAL=$(echo $LINE | awk '{print $1}')
	PID=$(echo $LINE | awk '{print $2}')
	OUTPUT=./RNA/mixcr/"${PID}"
fi

mixcr analyze shotgun -s hsa --starting-material "${TYPE}" "${LOC}"/"${NORMAL}"_1.fastq "${LOC}"/"${NORMAL}"_2.fastq "${OUTPUT}"
