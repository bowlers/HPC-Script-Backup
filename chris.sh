#!/bin/bash
#SBATCH --job-name=ChrisHRD
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/chris/chris-%A.err
#SBATCH --output=./logs/chris/chris-%A.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/maf

if [ -z "${1}" ]; then
        echo "Provide a study"
	exit 1
else
        STUDY="${1}"
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"

if [ ! -d "${DIR}" ]; then
        echo "Invalid study."
        exit 1
fi

cd "${DIR}"

Rscript /home/sbowler/biocore_lts/scott/chris/HRD.r

