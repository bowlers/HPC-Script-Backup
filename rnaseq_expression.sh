#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --partition=shared
##SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=18G
#SBATCH --error=./logs/rna/rnaseq-%A-%a.err
#SBATCH --output=./logs/rna/rnaseq-%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load lang/Perl

if [ -z "${1}" ]; then
	echo "You must provide a study number"
	exit 1
fi

STUDY="${1}"
if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
	echo "Could not locate study directory"
	exit 1
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA ]; then
	echo "Study directory does not contain RNA data"
	exit 1
fi

cd $HOME/biocore_lts/scott/"${STUDY}"/RNA/bam

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
	SLURM_ARRAY_TASK_ID=1
fi

FILE=../design.tsv
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
FNAME=$(echo $LINE | awk '{print $1}')

if [ -z "${FNAME}".bam ]; then
	echo "Could not locate fname"
	exit 1
fi

if [ ! -f "${FNAME}".bam.bai ]; then
	module load bio/SAMtools
	samtools index -b "${FNAME}".bam
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/stringtie ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/stringtie
fi

stringtie -p 2 -G $HOME/biocore_lts/scott/resources/GRCh38.gtf -e -B -o $HOME/biocore_lts/scott/"${STUDY}"/RNA/stringtie/"${FNAME}"/transcripts.gtf -A $HOME/biocore_lts/scott/"${STUDY}"/RNA/stringtie/"${FNAME}"/gene_abundances.tsv "${FNAME}".bam
