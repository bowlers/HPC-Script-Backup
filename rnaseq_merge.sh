#!/bin/bash
#SBATCH --job-name=RNAmerge
##SBATCH --partition=shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --error=./logs/rna/rnaseq-%A-%a.err
#SBATCH --output=./logs/rna/rnaseq-%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

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

cd $HOME/biocore_lts/scott/"${STUDY}"/RNA/htseq

if [ -f table.tsv ]; then
	rm table.tsv
fi

if [ -f final.tsv ]; then
	rm final.tsv
fi

header="GeneID "
x=0
for d in *.tsv ; do
        if [ "${d}" == "design.tsv" ]; then
                continue
        fi
        if [ "${x}" -eq 0 ]; then
                FIRST=${d}
        elif [ "${x}" -eq 1 ]; then
                join $FIRST $d > temp.tsv
        else
                join table.tsv $d > temp.tsv
        fi
	IFS='.' read -r -a nArray <<< "$d"
        iName=${nArray[0]}
        header+="${iName}"
        header+=" "
        if [ -f temp.tsv ]; then
                rm table.tsv
                mv temp.tsv table.tsv
        fi
        ((x++))
done

echo $header > header.txt
cat header.txt table.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > final.tsv
rm header.txt
rm table.tsv
