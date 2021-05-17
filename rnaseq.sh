#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --partition=shared
##SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --error=./logs/rna/rnaseq-%A-%a.err
#SBATCH --output=./logs/rna/rnaseq-%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/python3.7

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
	module load bio/SAMtools
	samtools index -b "${FNAME}".bam
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/htseq ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/htseq
fi

if [ ! -f $HOME/biocore_lts/scott/"${STUDY}"/RNA/htseq/"${FNAME}".tsv ]; then
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id "${FNAME}".bam $HOME/biocore_lts/scott/resources/hg38.gtf > ../htseq/"${FNAME}".tsv
fi
