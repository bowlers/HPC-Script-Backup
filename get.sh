#!/bin/bash

FILE=design.txt
touch compiled.tsv
echo -e "PID\tSample\ttotal_reads\tread_length\trepeat_threshold_set\trepeat_threshold_used\ttel_reads\tintratel_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\tTCAGGG_arbitrary_context_norm_by_intratel_reads\tTGAGGG_arbitrary_context_norm_by_intratel_reads\tTTGGGG_arbitrary_context_norm_by_intratel_reads\tTTCGGG_arbitrary_context_norm_by_intratel_reads\tTTTGGG_arbitrary_context_norm_by_intratel_reads\tATAGGG_arbitrary_context_norm_by_intratel_reads\tCATGGG_arbitrary_context_norm_by_intratel_reads\tCTAGGG_arbitrary_context_norm_by_intratel_reads\tGTAGGG_arbitrary_context_norm_by_intratel_reads\tTAAGGG_arbitrary_context_norm_by_intratel_reads\tTCAGGG_singletons_norm_by_all_reads\tTGAGGG_singletons_norm_by_all_reads\tTTGGGG_singletons_norm_by_all_reads\tTTCGGG_singletons_norm_by_all_reads\tTTTGGG_singletons_norm_by_all_reads\tATAGGG_singletons_norm_by_all_reads\tCATGGG_singletons_norm_by_all_reads\tCTAGGG_singletons_norm_by_all_reads\tGTAGGG_singletons_norm_by_all_reads\tTAAGGG_singletons_norm_by_all_reads" >> compiled.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET="./hunter/"
	TARGET+="${PID}"
	TARGET+="/"
	TARGET+="${PID}"_summary.tsv

	LINE1=$(sed -n '2p' "${TARGET}")
	LINE2=$(sed -n '3p' "${TARGET}")

	echo "${LINE1}" | tr ' ' '\t' >> compiled.tsv
	echo "${LINE2}" | tr ' ' '\t' >> compiled.tsv
done < "${FILE}"
