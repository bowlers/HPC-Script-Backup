#!/bin/bash 
#SBATCH --job-name=tmb-comb
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=36G
#SBATCH --cpus-per-task=2
#SBATCH --error=./logs/tmb/tmb-%A-%a.err
#SBATCH --output=./logs/tmb/tmb-%A-%a.out

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env

REFERENCE=/home/sbowler/biocore_lts/resources/GRCh38.fa

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
        echo "Invalid study."
        exit 1
fi

cd "${DIR}"

if [ ! -d ./snv ]; then
	mkdir snv
fi

if [ ! -f design.txt ]; then
	echo "Missing design.txt file containing TUMOR NORMAL PID format."
	exit 1
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" design.txt)
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

echo "${STUDY}"
echo "${PID}"


if [ -z "${2}" ]; then
	echo "Running Mutect2"
	gatk Mutect2 -R "${REFERENCE}" -I ./bam/"${TUMOR}".bam -I ./bam/"${NORMAL}".bam -normal "${NORMAL}" --germline-resource /home/sbowler/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz --panel-of-normals /home/sbowler/biocore_lts/scott/resources/1000g_pon.hg38.vcf.gz --f1r2-tar-gz ./snv/"${PID}".f1r2.tar.gz -O ./snv/"${PID}".somatic.vcf.gz
	sbatch --array="${SLURM_ARRAY_TASK_ID}" tmb.sh "${STUDY}" pileup
else
	ACTION="${2}"
	
	if [ "${ACTION}" == "pileup" ]; then
		echo "Running GetPileUpSummaries"
		gatk GetPileupSummaries -R "${REFERENCE}" -I ./bam/"${TUMOR}".bam -V /home/sbowler/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz -L /home/sbowler/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz -O ./snv/"${TUMOR}".pileups.table
		gatk GetPileupSummaries -R "${REFERENCE}" -I ./bam/"${NORMAL}".bam -V /home/sbowler/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz -L /home/sbowler/biocore_lts/scott/resources/af-only-gnomad.hg38.vcf.gz -O ./snv/"${NORMAL}".pileups.table
		sbatch --array="${SLURM_ARRAY_TASK_ID}" tmb.sh "${STUDY}" calc
	elif [ "${ACTION}" == "calc" ]; then
		echo "Running Calc Contamination"
		gatk CalculateContamination -I ./snv/"${TUMOR}".pileups.table -matched ./snv/"${NORMAL}".pileups.table -tumor-segmentaiton ./snv/"${PID}".segments.tsv -O ./snv/"${PID}".contamination.table
		sbatch --array="${SLURM_ARRAY_TASK_ID}" tmb.sh "${STUDY}" learn
	elif [ "${ACTION}" == "learn" ]; then
		echo "Running Learning"
		gatk LearnReadOrientationModel -I ./snv/"${PID}".f1r2.tar.gz -O ./snv/"${PID}".read-orientation-model.tar.gz
		sbatch --array="${SLURM_ARRAY_TASK_ID}" tmb.sh "${STUDY}" filter
	elif [ "${ACTION}" == "filter" ]; then
		echo "Running FilterMutect"
		gatk FilterMutectCalls -R "${REFERENCE}" -V ./snv/"${PID}".somatic.vcf.gz --contamination-table ./snv/"${PID}".contamination.table --tumor-segmentation ./snv/"${PID}".segments.tsv -O ./snv/"${PID}".filtered.vcf
		sbatch --array="${SLURM_ARRAY_TASK_ID}" tmb.sh "${STUDY}" funco
	elif [ "${ACTION}" == "funco" ]; then
		echo "Running Funcotator"
		gatk Funcotator -R "${REFERENCE}" -V ./snv/"${PID}".filtered.vcf -O ./snv/"${PID}".maf --output-file-format MAF --data-sources-path /home/sbowler/biocore_lts/scott/resources/funcotator_dataSources.v1.6.20190124s/ --ref-version hg38
	else
		echo "Invalid ACTION value."
		exit 1
	fi
fi
