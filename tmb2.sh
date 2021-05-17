#!/bin/bash 
#SBATCH --job-name=tmb2
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/tmb2-%A-%a.err
#SBATCH --output=./logs/tmb2-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${1}" ]; then
        echo "You must supply a study."
	exit 1
else
        STUDY="${1}"
fi

DIR="/home/sbowler/biocore_lts/scott/"
DIR+="${STUDY}"
DIR+="/"

if [ ! -d "${DIR}" ]; then
        echo "Invalid study.  Not located in /home/sbowler/biocore_lts/STUDY"
        exit 1
fi

if [ ! -d "${DIR}"/varscan ]; then
	echo "Study does not contain varscan data."
	exit 1
fi

cd "${DIR}"varscan
for i in *.snp.annot.mod.vcf; do
	bgzip $i
	tabix "${i}".gz
done
bcftools merge *.snp.annot.mod.vcf.gz --force-samples --merge all -O v -o Merged.snp.annot.mod.vcf
