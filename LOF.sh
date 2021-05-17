#!/bin/bash
#SBATCH --job-name=LOF
##SBATCH --partition=shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/LOF-%A-%a.err
#SBATCH --output=./logs/LOF-%A-%a.out

module load lang/Python

if [ -z "${1}" ]; then
	echo "You must supply a study."
	exit 1
else
	STUDY="${1}"
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
	echo "Invalid study"
	exit 1
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/varscan/lof/filtered/lof/tmp ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/varscan/lof/filtered/lof/tmp
fi

if [ -z "${2}" ]; then
	OUTPUT="final.tsv"
else
	OUTPUT="${2}"
fi

if [ ! -f $HOME/biocore_lts/scott/"${STUDY}"/varscan/lof/filtered/lof/unique.lst ]; then
	cd $HOME/biocore_lts/scott/"${STUDY}"/varscan/lof/filtered/lof
	cat *.txt | sort | uniq | cut -f1 | sort >> unique.lst
	cd $HOME/biocore_lts/scott/scripts
fi

python LOF.py -s "${STUDY}" -g unique.lst -o "${OUTPUT}"


