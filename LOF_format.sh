#!/bin/bash

module load lang/R
module load tools/nano

if [ -z "${1}" ]; then
	echo "You must provide a study"
	exit 1
else
	STUDY="${1}"
fi

cd $HOME/biocore_lts/scott/"${STUDY}"/varscan/lof/filtered

for i in `ls *.filtered.vcf`; do
	if [ ! -f ./lof/"${i}".LOF.txt ]; then
		Rscript format.R "${STUDY}" "${i}"
	fi
done

cd ./lof

for i in `ls *.LOF.txt`; do
	FILENAME=$( echo $i | awk -F. '{print $1}' )
	if [ ! -f "${FILENAME}".txt ]; then
		touch $FILENAME.txt
		while read LINE; do 
			NAME=$( echo $LINE | awk '{print $1}' )
			ID=$( echo $LINE | awk '{print $2}' )
			TRANSCRIPT=$( echo $LINE | awk '{print $3}' )
			FREQ=$( echo $LINE | awk '{print $4}' )

			TEMP="${ID}\t1"
			echo -e $TEMP >> $FILENAME.txt
		done < "${i}"
	fi
	sort -f "${FILENAME}".txt -o "${FILENAME}".sorted.txt
	rm "${FILENAME}".txt
	mv "${FILENAME}".sorted.txt "${FILENAME}".txt
	rm "${i}"
done
