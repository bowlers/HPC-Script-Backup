#!/bin/bash

while read p; do
	FILE=$(echo $p | awk '{print $1}')
	FILE2=$(echo $p | awk '{print $2}')
	if [ ! -f ./bam/"${FILE}".bam ]; then
		if [ ! -f ./bam/"${FILE2}".bam ]; then
			echo "$FILE and its partner does not exist"
		else
			echo "$FILE does not exist"
		fi
	fi

done < design.txt
