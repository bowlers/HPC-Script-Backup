#!/bin/bash

FILE=design.txt
touch msi_compiled.tsv
echo -e "PID\tTotal_sites\tNum_Somatic_sites\tMSI_score" >> msi_compiled.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET="./msisensor/"
	TARGET+="${PID}"
	TARGET+="/"
	TARGET+="output"

	AROW=$(sed -n '2p' "${TARGET}")

	TOTAL=$( echo "${AROW}" | awk '{print $1;}' )
	SOMATIC=$( echo "${AROW}" | awk '{print $2;}' )
	SCORE=$( echo "${AROW}" | awk '{print $3;}' )

	echo -e "$PID\t$TOTAL\t$SOMATIC\t$SCORE" >> msi_compiled.tsv
done < "${FILE}"
