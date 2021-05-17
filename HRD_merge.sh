#!/bin/bash
#SBATCH --job-name=HRDmerge
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/HRD/Merge-%A-%a.err
#SBATCH --output=./logs/HRD/Merge-%A-%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu


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

if [ -f HRD.tsv ]; then
	rm HRD.tsv
fi

touch HRD.tsv
# PID	LOH-HRD	Telomeric-AI	LST	HRD-Sum	cellularity	ploidy

echo -e "PID\tLOH-HRD\tTelomeric-AI\tLST\tHRD-Sum\tCellularity\tPloidy" >> HRD.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	HRD=./HRD/"${PID}".small.seqz._HRDresults.txt
	INFO=./HRD/"${PID}".small.seqz._info_seg.txt

	COUNTER=0
	LOH=0
	TELO=0
	LST=0
	SUM=0
	CELL=0
	PLO=0

	while read AROW; do
		((COUNTER+=1))
		if [ "${COUNTER}" -eq 1 ]; then
			continue
		fi
		LOH=$(echo $AROW | awk '{print $2}' )
		TELO=$(echo $AROW | awk '{print $3}' )
		LST=$(echo $AROW | awk '{print $4}' )
		SUM=$(echo $AROW | awk '{print $5}' )
	done < "${HRD}"

	COUNTER=0
	while read AROW; do
		((COUNTER+=1))
		if [ "${COUNTER}" -eq 1 ]; then
			continue
		fi
		CELL=$(echo $AROW | awk '{print $1}' )
		PLO=$(echo $AROW | awk '{print $2}' )
	done < "${INFO}"

	echo -e "$PID\t$LOH\t$TELO\t$LST\t$SUM\t$CELL\t$PLO" >> HRD.tsv
done < "${FILE}"
