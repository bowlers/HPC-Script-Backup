#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --partition=shared
##SBATCH --partition=kill-exclusive
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --error=./logs/rnaseq-%A-%a.err
#SBATCH --output=./logs/rnaseq-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

# ---------------------------------------------------------------------------------------------------------
#
# Usage: Called from RNA_batch.sh
# Version: 1.02
# Date: July 23, 2020
# Author: Scott Bowler
#
# ---------------------------------------------------------------------------------------------------------

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/env
module load lang/Perl
module load bio/SAMtools

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

AROW=1
for i in `ls *.bam | sed -n 's/\.bam$//p'`; do
	if [ "${SLURM_ARRAY_TASK_ID}" -eq "${AROW}" ]; then
		FNAME="${i}"
		break
	else
		((AROW++))
	fi
done

if [ -z "${FNAME}".bam ]; then
	echo "Could not locate fname"
	exit 1
fi

if [ ! -f "${FNAME}".bam.bai ]; then
	samtools index -@ 8 "${FNAME}".bam
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete/"${FNAME}" ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete/"${FNAME}"
fi

stringtie -p 8 -G $HOME/biocore_lts/scott/resources/hg38.ensGene.gtf -e -B -o ../complete/"${FNAME}"/transcripts.gtf -A ../complete/"${FNAME}"/gene_abundances.tsv "${FNAME}".bam

