#!/bin/bash
#SBATCH --job-name=CNV
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --error=./logs/cnv/cnv-%A-%a.err
#SBATCH --output=./logs/cnv/cnv-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
module load bio/SAMtools

RES=$HOME/biocore_lts/scott/resources/GRCh38.fa
REF=$HOME/biocore_lts/scott/resources

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

# Locate design file so that we can pull normal and tumor SRRs and PID
if [ ! -f "${STUDY_DIR}"/design.txt ]; then
        echo No design file within study. >&2
        exit 1
fi

if [ ! -d "${STUDY_DIR}"/loh ]; then 
	mkdir "${STUDY_DIR}"/loh
fi


# Design exists, so lets pull out the tumor and normal SRRs and PID
DESIGN="${STUDY_DIR}/design.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${DESIGN}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

cd "${STUDY_DIR}"/loh

if [ ! -f "${STUDY_DIR}"/loh/"${PID}".copynumber ]; then
	source activate /home/sbowler/biocore_lts/scott/envir/env
	samtools mpileup -q 1 -f "${RES}" "${STUDY_DIR}"/bam/"${NORMAL}".bam "${STUDY_DIR}"/bam/"${TUMOR}".bam > "${PID}".pileup
	varscan copynumber "${PID}".pileup "${PID}" --mpileup 1
	varscan copyCaller "${PID}".copynumber --output-file "${PID}".copynumber.called 
fi

if [ ! -f "${STUDY_DIR}"/loh/"${PID}".copynumber.called.output ]; then
	source activate /home/sbowler/biocore_lts/scott/envir/R
	Rscript /home/sbowler/biocore_lts/scott/scripts/CNV.R "${STUDY}" "${PID}"
fi

if [ ! -f "${STUDY_DIR}"/loh/"${PID}".summary.txt ]; then
	module load lang/Perl
        source activate /home/sbowler/biocore_lts/scott/envir/R
	tail -n +2 "${PID}".copynumber.called.output > "${PID}".noheader
	perl /home/sbowler/biocore_lts/scott/scripts/mergeSegments.pl "${PID}".noheader --ref-arm-sizes "${REF}"/hg38.refarm.txt --output-basename "${PID}"
fi
