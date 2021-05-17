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

if [ ! -d "${DIR}"/varscan ]; then
	echo "Study does not contain varscan data."
	exit 1
fi

cd "${DIR}"varscan
AROW=1
for i in `ls *.snp.vcf.gz | sed -n 's/\.snp.vcf.gz$//p'`;
        do
                if [ "${AROW}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
                        FNAME="${i}"
                        break
                fi
                ((AROW++))
        done

# bgzip "${FNAME}".snp.vcf
# tabix "${FNAME}".snp.vcf.gz
# bcftools norm -m-both -o "${FNAME}".snp.step1.vcf "${FNAME}".snp.vcf.gz
#  bcftools norm -f $HOME/biocore_lts/scott/resources/GRCh38.fa -o "${FNAME}".snp.step2.vcf "${FNAME}".snp.step1.vcf
#  $HOME/biocore_lts/scott/annovar/convert2annovar.pl -format vcf4 "${FNAME}".snp.step2.vcf > "${FNAME}".snp.avinput
# cd $HOME/biocore_lts/scott/annovar
# ./table_annovar.pl "${DIR}"varscan/"${FNAME}".snp.avinput humandb/ -buildver hg38 -out "${DIR}"varscan/"${FNAME}".myanno -remove -protocol refGene -operation g -nastring . -vcfinput

cp header.txt "${FNAME}".snp.annot.mod.vcf
bcftools annotate -a $HOME/biocore_lts/scott/resources/hg38_RefSeq.sorted.bed.gz -c CHROM,FROM,TO,GENE -h $HOME/biocore_lts/scott/scripts/header.txt "${FNAME}".snp.vcf.gz > "${FNAME}".snp.annot.vcf
bgzip "${FNAME}".snp.annot.vcf
tabix "${FNAME}".snp.annot.vcf.gz
bcftools query -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23 -f '%CHROM\t%POS\t%REF\t%ALT\t%SS\t%GENE[\t%SAMPLE=%GT]\n' "${FNAME}".snp.annot.vcf.gz >> "${FNAME}".snp.annot.mod.vcf
