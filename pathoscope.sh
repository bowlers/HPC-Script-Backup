#!/bin/bash
#SBATCH --job-name=patho
#SBATCH --partition=kill-exclusive
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=19
#SBATCH --mem=0
#SBATCH --error=./logs/patho/patho-$1-%A-%a.err
#SBATCH --output=./logs/patho/patho-$1-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/patho
module load bio/Bowtie2
module load lang/Python/2.7.15-intel-2018.5.274

if [ -z "${1}" ]; then
	echo Must provide a study. >&2
	exit 1
else
	STUDY="${1}"
	FILE_DIR="/home/sbowler/biocore_lts/scott/"
	FILE_DIR+="${STUDY}"
	FILE_DIR+="/"

	if [ ! -d "${FILE_DIR}" ]; then
		echo "$1 does not exist at $FILE_DIR."
		exit 1
	fi
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
	SLURM_ARRAY_TASK_ID=1
fi

PREFIX="/home/sbowler/biocore_lts/scott/"
PREFIX+="${STUDY}"
PREFIX+="/fastq/"
AROW=1
for i in `ls "${PREFIX}"*_1.fastq | sed -n 's/\_1.fastq$//p'`;
        do
                if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
                        i=${i#"$PREFIX"}
                        FNAME="${i}"
                        break
                fi
                ((AROW++))
        done

COMPLETE="${FILE_DIR}"
COMPLETE+="microbiome"
FILE_DIR+="fastq/"

echo "${STUDY}"
echo "${FNAME}"

if [ ! -f "${COMPLETE}"/"${FNAME}"-sam-report.tsv ]; then
	#create directory within DATA that includes study for longterm sorting between projects
	if [ ! -d "${COMPLETE}" ]; then
		mkdir "${COMPLETE}"
	fi
	if [ ! -d "${COMPLETE}"/"${FNAME}" ]; then
		mkdir "${COMPLETE}"/"${FNAME}"
	fi

	cd /home/sbowler/biocore_lts/scott/microbiome/references
	pathoscope MAP -numThreads 19 -1 "${FILE_DIR}${FNAME}"_1.fastq -2 "${FILE_DIR}${FNAME}"_2.fastq -targetIndexPrefixes Bacteria_0,Bacteria_1,Bacteria_2,Bacteria_3,Bacteria_4,Bacteria_5,Bacteria_6,Bacteria_7,Bacteria_8,Bacteria_9,Bacteria_10,Bacteria_11,Bacteria_12,Bacteria_13,Bacteria_14,Bacteria_15,Bacteria_16 -filterIndexPrefixes hg38,Phix174 -outDir "${COMPLETE}"/"${FNAME}" -outAlign "${FNAME}"_cpu.sam -expTag "${FNAME}"
	pathoscope ID -alignFile "${COMPLETE}"/"${FNAME}"/"${FNAME}"_cpu.sam -fileType sam -outDir "${COMPLETE}"/"${FNAME}" -expTag "${FNAME}" -thetaPrior 1000000
	mv "${COMPLETE}"/"${FNAME}"/"${FNAME}"-sam-report.tsv "${COMPLETE}"
	rm -r "${COMPLETE}"/"${FNAME}"/
fi
