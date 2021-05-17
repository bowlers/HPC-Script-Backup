#!/bin/bash
#SBATCH --job-name=DESeq2
##SBATCH --partition=shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --error=./logs/deseq/DESeq2r-%A-%a.err
#SBATCH --output=./logs/deseq/DESeq2r-%A-%a.out

module load lang/R

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

STUDY="${1}"

Rscript DESeq2.R "${STUDY}"
