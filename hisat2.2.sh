#!/bin/bash 
#SBATCH --job-name=hisat2 
##SBATCH --partition=shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=72G
#SBATCH --error=./logs/hisat2/hisat2-%A_%a.err
#SBATCH --output=./logs/hisat2/hisat2-%A_%a.out
#
#-----------------------------------------------------------------------------------------
#  Usage
#  sbatch --array=1-[max number of samples] hisat2 [option]
#-----------------------------------------------------------------------------------------
#  Local Definitions

REF_DIR=$HOME/biocore_lts/scott/resources

#-----------------------------------------------------------------------------------------

module load bio/HISAT2
module load lang/Python
module load bio/SAMtools

if [ -z "${1}" ]; then
	echo "You must supply a study name"
	exit 1
fi

STUDY="${1}"

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
	echo "Could not locate study directory"
	exit 1
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/fastq ]; then
	echo "Study directory does not seem to contain RNA files."
	exit 1
fi

cd $HOME/biocore_lts/scott/"${STUDY}"/RNA/fastq

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/summary ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/summary
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/bam ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/RNA/bam
fi

AROW=1
if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
	SLURM_ARRAY_TASK_ID=1
fi

for i in `ls *_1.fastq | sed -n 's/\_1.fastq$//p'`; do
	if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
		FNAME="${i}"
		break
	else
		((AROW++))
	fi
done

	hisat2 -p 10 --dta -x "${REF_DIR}"/GRCh38 -1 "${FNAME}"_1.fastq -2 "${FNAME}"_2.fastq -S "${FNAME}".sam --summary-file ../summary/"${FNAME}".txt
	samtools sort -@ 10 -o ../bam/"${FNAME}".bam "${FNAME}".sam
	samtools index -@ 10 "${FNAME}".bam
#	rm "${FNAME}".sam
