#!/bin/bash

for FILE in `ls *.*`; do
	if [ "${FILE}" == "edit_filenames.sh" ]; then
		continue
	fi
	EXT=${FILE: -8}

	if [ "${EXT}" != ".bam.bai" ]; then
		EXT=${FILE: -4}
	fi
	NEW_FILE=${FILE#"wo29074_ngs_dna_wes_sureselect_"}
	NEW_FILE=${NEW_FILE#"WO29074_ngs_dna_wes_sureselect_"}
	NEW_FILE=${NEW_FILE%_*}
	NEW_FILE+=$EXT

	mv "${FILE}" "${NEW_FILE}"
done
