#!/bin/bash
#SBATCH --job-name=tRF
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=54G
#SBATCH --error=./logs/tRF-%A-%a.err
#SBATCH --output=./logs/tRF-%A-%a.out
#SBATCH --mail-type=
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Python/3.7.2-intel-2018.5.274

if [ -z "${1}" ]; then
        echo Must provide a study. >&2
        exit 1
fi
STUDY="${1}"

python RanFor.py "${STUDY}"


