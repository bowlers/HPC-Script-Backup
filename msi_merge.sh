#!/bin/bash

FILE=design.txt

if [ -z "${1}" ]; then
	echo Provide a study number. >&2
	exit 1
fi

STUDY="${1}"

if [ ! -d /home/sbowler/biocore_lts/scott/"${STUDY}" ]; then
	echo Could not locate study directory. >&2
	exit 1
fi

cd /home/sbowler/biocore_lts/scott/"${STUDY}"

touch msi.tsv
# PID 	Total_Number_of_Sites   Number_of_Somatic_Sites %
echo -e "PID\tTotal_Number_of_Sites\tNumber_of_Somatic_Sites\t%" >> msi.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET="./msisensor/"
	TARGET+="${PID}"
	TARGET+="/output"

	COUNTER=1
	while read AROW; do
		if [ "${COUNTER}" -gt 1 ]; then
			echo -e "${PID}\t${AROW}" >> msi.tsv
		fi
		((COUNTER+=1))
	done < "${TARGET}"
done < "${FILE}"
