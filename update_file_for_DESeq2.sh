#!/bin/bash

while read -r line; do
	FNAME=$(echo $line | awk '{print $1}')
	STATUS=$(echo $line | awk '{print $2}')

	VAR="non"
	if [ $STATUS -eq 0 ]; then
		VAR="resp"
	fi

	mv ./htseq/"${FNAME}".tsv ./htseq/"${FNAME}"."${VAR}".tsv
done < design.tsv
