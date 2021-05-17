#!/bin/bash
#SBATCH --job-name=LOFilter
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --error=./logs/vscan/LOFilter-%A-%a.err
#SBATCH --output=./logs/vscan/LOFilter-%A-%a.out

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

STUDY="${1}"
if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
        echo "Could not locate study directory"
        exit 1
fi
STUDY_DIR=$HOME/biocore_lts/scott/"${STUDY}"

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

if [ ! -d "${STUDY_DIR}"/varscan ]; then
        echo "Could not locate somatic mutation data."
	exit 1
fi

cd "${STUDY_DIR}"/varscan

if [ ! -d ./lof ]; then
	echo "Rerun varscan.sh, could not locate LOF snpeff data."
	exit 1
fi
cd ./lof

if [ ! -d ./filtered ]; then
	mkdir filtered
fi

for f in `ls *.vcf`; do
	bgzip -c "${f}" > "${f}".gz
	tabix "${f}".gz
	bcftools filter -i 'INFO/LOF!="."' "${f}".gz > ./filtered/"${f}".filtered.vcf
	rm "${f}".gz
	rm "${f}".gz.tbi
done
