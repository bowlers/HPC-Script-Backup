#!/bin/bash
#SBATCH --job-name=NBS
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --partition=kill-shared
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --error=./logs/NBS-%A-%a.err
#SBATCH --output=./logs/NBS-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/nbs
export PATH=$HOME/biocore_lts/scott/pyNBS:$PATH

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
	echo "Please supply a study name"
	exit 1
fi

if [ -z "${2}" ]; then
	echo "No k value supplied as variable input 2, defaulting to 2."
	K=2
else
	K="${2}"
fi

STUDY="${1}"
if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
	echo "Could not locate study directory."
	exit 1
fi

DIR=$HOME/biocore_lts/scott/"${STUDY}"
if [ ! -d "${DIR}"/NBS ]; then
	echo "Study does not have NBS directory of data."
	exit 1
fi

if [ ! -f "${DIR}"/NBS/network.tsv ]; then
	echo "Study NBS directory does not contain gene network (network.tsv)."
	exit 1
fi
NETWORK="${DIR}"/NBS/network.tsv

if [ ! -f "${DIR}"/NBS/matrix.tsv ]; then
	echo "Study NBS directory does not contain mutation matrix (matrix.tsv)."
	exit 1
fi
MATRIX="${DIR}"/NBS/matrix.tsv

if [ ! -d "${DIR}"/NBS/output ]; then
	mkdir "${DIR}"/NBS/output
fi

python ../pyNBS/run_pyNBS.py -params ../pyNBS/params.csv -k "${K}" -o "${DIR}"/NBS/output/ "${MATRIX}" "${NETWORK}"
