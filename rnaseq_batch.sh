#!/bin/bash

if [ -z "$1" ]
then
	echo "Syntax: sh rnaseq_batch [number] [option]"
	echo "Options: -nodownload"
else
	COUNT=$1
	NODOWNLOAD=false
	if [[ -n "${2}" && "${2}" = "-nodownload" ]]; then
		NODOWNLOAD=true
	fi

	if "${NODOWNLOAD}" = "true"; then
		jid1=$(sbatch --array=1-"${COUNT}" rnaseq_expression.sh -nodownload | cut -f 4 -d' ')
	else
		jid1=$(sbatch --array=1-"${COUNT}" rnaseq_expression.sh | cut -f 4 -d' ')
	fi

	jid2=$(sbatch --dependency=afterok:"$jid1" rnaseq_complete.sh)
fi
