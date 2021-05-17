#!/bin/bash
#SBATCH --job-name=CNVr
#SBATCH --partition=shared
##SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --error=./logs/CNVR-%A-%a.err
#SBATCH --output=./logs/CNVR-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/env
module load bio/SAMtools
module load lang/R/3.5.1-intel-2018.5.274-Python-2.7.15
export R_LIBS=/home/sbowler/biocore_lts/scott/Rlib


RES=$HOME/biocore_lts/scott/resources/hg38.fa

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

STUDY="${1}"

Rscript CNV.R "${STUDY}"
