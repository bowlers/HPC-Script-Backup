#!/bin/bash 
#SBATCH --job-name=annotate
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/tmb-anno-%A-%a.err
#SBATCH --output=./logs/tmb-anno-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/env

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
        echo "Invalid study.  Not located in /home/sbowler/biocore_lts/STUDY"
        exit 1
fi

if [ ! -d "${DIR}"/varscan ]; then echo "Study does not contain varscan data." exit 1
fi

cd "${DIR}"varscan
if [ ! -f "${DIR}"varscan/Merged.snp.vcf.gz ]; then
	bgzip Merged.snp.vcf
	tabix Merged.snp.vcf.gz
fi

# bgzip Merged.snp.vcf
bcftools annotate -a $HOME/biocore_lts/scott/resources/refGene.hg38.sorted.bed -c CHROM,FROM,TO,GENE -h $HOME/biocore_lts/scott/scripts/header.txt Merged.snp.vcf.gz > Merged.snp.annot.vcf
# bgzip Merged.snp.annot.vcf
# tabix Merged.snp.annot.vcf.gz
# bcftools query -H -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23 -f '%CHROM\t%POS\t%REF\t%ALT\t%SS\t%GENE[\t%SAMPLE=%GT]\n' Merged.snp.annot.vcf.gz > Merged.snp.annot.mod.vcf
