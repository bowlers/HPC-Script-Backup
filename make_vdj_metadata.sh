#!/bin/bash
if [ -f metadata.txt ]; then
	rm metadata.txt
	touch metadata.txt
fi

echo -e "#file.name\tsample.id\tfilter\tresponse" > metadata.txt

for i in `ls *.ALL.txt | sed -n 's/\.clonotypes.ALL.txt$//p'`; do
	TEXT=""
	TEXT+=$i
	TEXT+=".txt\t"
	TEXT+=$i
	TEXT+="\t0\t"
	# Add response values
	TEXT+="NR"
	echo -e $TEXT >> metadata.txt
done
	
