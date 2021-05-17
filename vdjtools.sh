#!/bin/bash

module load lang/Java
VDJTOOLS="java -Xmx20G -jar /home/sbowler/biocore_lts/scott/vdjtools/vdjtools.jar"

if [ -z "${1}" ]; then
	echo "Provide a study"
	exit 1
else
	STUDY="${1}"
fi

if [ ! -d /home/sbowler/biocore_lts/scott/"${STUDY}" ]; then
	echo "Invalid study directory"
	exit 1
else
	cd /home/sbowler/biocore_lts/scott/"${STUDY}" ];
fi

if [ ! -d ./mixcr ]; then
	echo "Could not locate mixcr data."
	exit 1
else
	cd ./mixcr
fi

if [ ! -f ./metadata.txt ]; then
	echo "Could not locate metadata file."
	exit 1
elif [ ! -d ./converted ]; then
	echo "Could not located converted data."
	exit 1
else
	cd ./converted
fi

if [ ! -d ./out ]; then
	mkdir ./out
fi

$VDJTOOLS CalcBasicStats -m ../metadata.txt ./out
$VDJTOOLS CalcSpectratype -m ../metadata.txt ./out

# -p for plotting, -f specifies metadata column for coloring, -n tells that factor is continous (ours is categorical)
mkdir ./out/segUsage
$VDJTOOLS CalcSegmentUsage -m ../metadata.txt -p -f response ./out/segUsage

while read LINE; do
	FILE=$(echo $LINE | awk '{print $1}')
	PID=$(echo $LINE | awk '{print $2}')
	mkdir ./out/$PID
	$VDJTOOLS PlotFancySpectratype $FILE ./out/$PID
	$VDJTOOLS PlotSpectratypeV $FILE ./out/$PID
	$VDJTOOLS PlotFancyVJUsage $FILE ./out/$PID
done < ../metadata.txt

mkdir ./out/dist
# Characterize similiarties between repertoires
$VDJTOOLS CalcPairwiseDistances -m ../metadata.txt ./out/dist

# Sample Overlapping
# Uncomment if you can justify comparing T v N.  Does WXS tumor make sense?  What about using T RNAseq data and N WXS?
# Maybe try this in the Native Hawaiian cohort since they have T/N matched RNA?
# Discuss w Ved about feasibility
# Use design file to obtain T and N filenames (note T samples have NOT been ran through mixCR as of 2021.03.09) SAB
# ---------------------------------------------------------------------------------------------------------------------
# NOTE: this will not successfully run until the samples are rerun with mixcr to change the names to their SRR 
#       designations.  Currently they are coded as PID since only N was done.  Minor change, but DONT FORGET TO DO IT!
# ---------------------------------------------------------------------------------------------------------------------
#
#
# while read LINE; do
#	TUMOR=$(echo $LINE | awk '{print $1}')
#	NORMAL=$(echo $LINE | awk '{print $2}')
#	PID=$(echo $LINE | awk '{print $3}')
#
#	mkdir ./out/$PID/overlap
#	Overlap two replicate samples from same donor
#	$VDJTOOLS OverlapPair -p $TUMOR.txt $NORMAL.txt ./out/$PID/overlap
# done < ../../design.txt

# Sample clustering
$VDJTOOLS ClusterSamples -p -f response -1 sample.id ./out/dist ./out/dist.response
# Jensen-Shannon divergence 
$VDJTOOLS ClusterSamples -p -e vJSD -f age -n -1 sample.id ./out/dist ./out/dist.age

# Remove contamination
mkdir ./out/decon
$VDJTOOLS Decontaminate -m ../metadata.txt -c ./out/decon

# Down-sample datasets to 10,000 reads
mkdir ./out/down
$VDJTOOLS Downsample -m ../metadata.txt -c -x 10000 ./out/down

# Filter out non-cloning types
mkdir ./out/filtered
$VDJTOOLS FilterNonFunctional -m ../metadata.txt -c ./out/filtered

# Join samples into a single matrix
mkdir ./out/joined
$VDJTOOLS JoinSamples -p -m ../metadata.txt ./out/joined

# Everyone likes a Pool
mkdir ./out/pooled
$VDJTOOLS PoolSamples -m ../metadata.txt ./out/pooled

# Annotate
mkdir ./out/annot
$VDJTOOLS Annotate -m ../metadata.txt ./out/annot
