#!/bin/bash
#SBATCH --job-name=rainfall
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/rainfall-%A-%a.err
#SBATCH --output=./logs/rainfall-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/maf

echo "${STUDY}"
echo "${PID}"

Rscript /home/sbowler/biocore_lts/scott/scripts/rainfall.R
