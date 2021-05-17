#!/bin/bash 
#SBATCH --job-name=varscan
#SBATCH --partition=kill-shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --error=./logs/TMB/vTMB-%A-%a.err
#SBATCH --output=./logs/TMB/vTMB-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/env
module load bio/SAMtools
module load lang/Perl
module load bio/BioPerl

RES=$HOME/biocore_lts/scott/U54/Homo_sapiens.GRCh38.chr.toplevel.fa

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
	mkdir "${STUDY_DIR}"/varscan
fi

cd "${STUDY_DIR}"
FILE="${DIR}"
FILE+=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" "${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')
FNAME="${PID}"

if [ ! -d "${STUDY_DIR}"/varscan/completed ]; then
        mkdir "${STUDY_DIR}"/varscan/completed
fi
if [ ! -d "${STUDY_DIR}"/varscan/filtered ]; then
        mkdir "${STUDY_DIR}"/varscan/filtered
fi
if [ ! -d "${STUDY_DIR}"/varscan/maf ]; then
        mkdir "${STUDY_DIR}"/varscan/maf
fi

cd "${STUDY_DIR}"/varscan

echo "${FNAME}"

if [ ! -f ./completed/"${FNAME}".snp.vcf ]; then
	samtools mpileup -q 1 -f "${RES}" ../bam/"${NORMAL}".bam ../bam/"${TUMOR}".bam > "${PID}".pileup
	varscan somatic "${STUDY_DIR}"/varscan/"${PID}".pileup "${STUDY_DIR}"/varscan/"${PID}" --somatic-p-value 0.05 --p-value 0.99 --output-vcf 1 --mpileup 1
	mv "${FNAME}".snp.vcf ./completed
	mv "${FNAME}".indel.vcf ./completed
fi

if [ ! -f ./maf/"${FNAME}".vep.maf ]; then
	varscan somaticFilter ./completed/"${PID}".snp.vcf --indel-file ./completed/"${PID}".indel.vcf --output-file ./filtered/"${PID}".filtered.snp.vcf
	perl $HOME/biocore_lts/scott/vcf2maf/vcf2maf.pl --input-vcf ./filtered/"${PID}".filtered.snp.vcf --output-maf ./maf/"${PID}".vep.maf --tumor-id "${TUMOR}" --normal-id "${NORMAL}" --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --ref-fasta "${RES}" --vep-path $HOME/biocore_lts/scott/envir/env/bin --vep-data $HOME/biocore_lts/scott/vep --ncbi-build GRCh38
fi

