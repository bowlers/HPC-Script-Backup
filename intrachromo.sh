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

touch intra_tumor.tsv
touch intra_control.tsv
# chr     band    reads_with_pattern      TTAGGG  CCCTAA  TGAGGG  CCCTCA  TCAGGG  CCCTGA  TTGGGG  CCCCAA  other

echo -e "PID\tchr\tband\treads_with_pattern\tTTAGGG\tCCCTAA\tTGAGGG\tCCCTCA\tTCAGGG\tCCCTGA\tTTGGGG\tCCCCAA\tother" >> intra_tumor.tsv
echo -e "PID\tchr\tband\treads_with_pattern\tTTAGGG\tCCCTAA\tTGAGGG\tCCCTCA\tTCAGGG\tCCCTGA\tTTGGGG\tCCCCAA\tother" >> intra_control.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET="./hunter/"
	TARGET+="${PID}"
	TARGET2="${TARGET}"
	TARGET+="/tumor_TelomerCnt_"
	TARGET+="${PID}"
	TARGET+=/"${PID}"_spectrum.tsv
	TARGET2+=/control_TelomerCnt_"${PID}"/"${PID}"_spectrum.tsv

	COUNTER=1
	while read AROW; do
		if [ "${COUNTER}" -gt 1 ]; then
			echo -e "${PID}\t${AROW}" >> intra_tumor.tsv
		fi
		((COUNTER+=1))
	done < "${TARGET}"

	COUNTER=1
	while read AROW; do
		if [ "${COUNTER}" -gt 1 ]; then
			echo -e "${PID}\t${AROW}" >> intra_control.tsv
		fi
		((COUNTER+=1))
	done < "${TARGET2}"
done < "${FILE}"
