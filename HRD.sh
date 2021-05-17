#!/bin/bash
#SBATCH --job-name=HRD
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=3
#SBATCH --error=./logs/HRD/HRD-%A-%a.err
#SBATCH --output=./logs/HRD/HRD-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env
module load lang/R/3.5.1-intel-2018.5.274-Python-2.7.15

#R library definitions
if [ -n $R_LIBS ]; then
        export R_LIBS=/home/sbowler/biocore_lts/scott/Rlib:$R_LIBS
else
        export R_LIBS=/home/sbowler/biocore_lts/scott/R_CNV
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
        SLURM_ARRAY_TASK_ID=1
fi

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

if [ ! -d "${DIR}"/HRD ]; then
	mkdir HRD
fi

FILE=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

cd HRD

echo "${STUDY}"
echo "${PID}"

if [ ! -e "${PID}".small.seqz._HRDresults.txt ] ; then
	if [ ! -e "${PID}".small.seqz.gz ] ; then
		sequenza-utils bam2seqz -n ../bam/"${NORMAL}".bam -t ../bam/"${TUMOR}".bam --fasta ../../resources/GRCh38.fa -gc ../../resources/GRCh38.gc50.wig.gz -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o "${PID}".out.seqz.gz
		sequenza-utils seqz_binning --seqz "${PID}".out.seqz.gz -w 50 -o "${PID}".small.seqz.gz
	#	rm "${PID}".out.seqz.gz
	#	rm "${PID}".out.seqz.gz.tbi
	fi
	Rscript /home/sbowler/biocore_lts/scott/scripts/HRD.R "${STUDY}" "${PID}"
fi
