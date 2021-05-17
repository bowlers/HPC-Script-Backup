#!/bin/bash
#SBATCH --job-name=RNAseq
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --error=./logs/rComp-%A_%a.err
#SBATCH --output=./logs/rComp-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load lang/Anaconda3
module load lang/Perl
source activate /home/sbowler/biocore_lts/scott/env/

if [ -z "${1}" ]; then
        echo "You must provide a study number"
        exit 1
fi

STUDY="${1}"
if [ ! -d $HOME/biocore_lts/scott/"${STUDY}" ]; then
        echo "Could not locate study directory"
        exit 1
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete ]; then
	echo "Could not locate completed RNA directory."
	exit 1
fi

cd $HOME/biocore_lts/scott/"${STUDY}"/RNA/complete

list=""
x=0

for d in */ ; do
	if [[ "$x" -eq 1 ]]; then
		list+=,
	fi
	iName="${d%?}"
	list+=$iName
	x=1
done

echo $list

$HOME/biocore_lts/scott/scripts/stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs="$list" --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

