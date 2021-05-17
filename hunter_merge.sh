#!/bin/bash
#SBATCH --job-name=HuntMerg
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

if [ -f TLS.tsv ]; then
        rm TLS.tsv
fi

touch TLS.tsv

FILE=design.txt
echo -e "PID\tSample\ttotal_reads\tread_length\trepeat_threshold_set\trepeat_threshold_used\ttel_reads\tintratel_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\tTCAGGG_arbitrary_context_norm_by_intratel_reads\tTGAGGG_arbitrary_context_norm_by_intratel_reads\tTTGGGG_arbitrary_context_norm_by_intratel_reads\tTTCGGG_arbitrary_context_norm_by_intratel_reads\tTTTGGG_arbitrary_context_norm_by_intratel_reads\tATAGGG_arbitrary_context_norm_by_intratel_reads\tCATGGG_arbitrary_context_norm_by_intratel_reads\tCTAGGG_arbitrary_context_norm_by_intratel_reads\tGTAGGG_arbitrary_context_norm_by_intratel_reads\tTAAGGG_arbitrary_context_norm_by_intratel_reads\tTCAGGG_singletons_norm_by_all_reads\tTGAGGG_singletons_norm_by_all_reads\tTTGGGG_singletons_norm_by_all_reads\tTTCGGG_singletons_norm_by_all_reads\tTTTGGG_singletons_norm_by_all_reads\tATAGGG_singletons_norm_by_all_reads\tCATGGG_singletons_norm_by_all_reads\tCTAGGG_singletons_norm_by_all_reads\tGTAGGG_singletons_norm_by_all_reads\tTAAGGG_singletons_norm_by_all_reads" >> TLS.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET="./hunter/"
	TARGET+="${PID}"
	TARGET+="/"
	TARGET+="${PID}"_summary.tsv

	LINE1=$(sed -n '2p' "${TARGET}")
	LINE2=$(sed -n '3p' "${TARGET}")

	echo "${LINE1}" | tr ' ' '\t' >> TLS.tsv
	echo "${LINE2}" | tr ' ' '\t' >> TLS.tsv
done < "${FILE}"
