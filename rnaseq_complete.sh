#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --error=./logs/rComp-%A_%a.err
#SBATCH --output=./logs/rComp-%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
source activate /home/sbowler/biocore_lts/scott/envir/python3.7
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

cd $HOME/biocore_lts/scott/"${STUDY}"/RNA/

file=design.tsv
list=""
x=0
while read -r line; do
	if [[ "$x" -eq 1 ]]; then
		list+=,
	fi
	iName=$( echo $line | awk '{print $1}' )
	list+=$iName
	x=1
done < $file

echo $list

cd ./stringtie

$HOME/biocore_lts/scott/scripts/stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs="$list" --transcript_matrix_file=transcript_tpms_all_samples.tsv --gene_matrix_file=gene_tpms_all_samples.tsv
