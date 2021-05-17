#!/bin/bash
#SBATCH --job-name=LOHmerge
##SBATCH --partition=kill-shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=03-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --error=./logs/loh/loh-%A-%a.err
#SBATCH --output=./logs/loh/loh-%A-%a.out
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

if [ -f ./loh.tsv ]; then
	rm loh.tsv
fi

touch loh.tsv
# PID	large-scale	focal	large-scale	focal	amp_regions	del_regions
echo -e "PID\tLarge-Scale-amp\tFocal\tLarge-Scale-del\tFocal\tAmp_Regions\tDel_Regions" >> loh.tsv

while read LINE; do
	PID=$(echo $LINE | awk '{print $3}')
	TARGET=./loh/"${PID}".events.tsv

	COUNTER=0
	LARGE_AMP=0
	FOCAL_AMP=0
	LARGE_DEL=0
	FOCAL_DEL=0
	AMP_REG=""
	DEL_REG=""

	while read AROW; do
		((COUNTER+=1))
		if [ "${COUNTER}" -eq 1 ]; then
			continue
		fi

		CHROM=$(echo $AROW | awk '{print $1}' )
                START=$(echo $AROW | awk '{print $2}' )
                STOP=$(echo $AROW | awk '{print $3}' )
                SEG_MEAN=$(echo $AROW | awk '{print $4}' )
                NUM_SEGS=$(echo $AROW | awk '{print $5}' )
                NUM_MARK=$(echo $AROW | awk '{print $6}' )
                P_VALUE=$(echo $AROW | awk '{print $7}' )
                EVENT=$(echo $AROW | awk '{print $8}' )
                E_SIZE=$(echo $AROW | awk '{print $9}' )
                SIZE_CLASS=$(echo $AROW | awk '{print $10}' )
                CHROM_ARM=$(echo $AROW | awk '{print $11}' )
                ARM_FRAC=$(echo $AROW | awk '{print $12}' | rev | cut -c 2- | rev )
                CHROM_FRAC=$(echo $AROW | awk '{print $13}' | rev | cut -c 2- | rev )
		ARM_FRAC=${ARM_FRAC%.*}
		CHROM_FRAC=${CHROM_FRAC%.*}

		if [[ $ARM_FRAC -ge 90 || $CHROM_FRAC -ge 90 ]]; then
			continue
		fi

		if [[ $EVENT == "amplification" ]]; then
			if [[ $SIZE_CLASS == "large-scale" ]]; then
				((LARGE_AMP+=1))
				AMP_REG+=$CHROM_ARM
				AMP_REG+=","
			else
				((FOCAL_AMP+=1))
			fi
		else
                        if [[ $SIZE_CLASS == "large-scale" ]]; then
                                ((LARGE_DEL+=1))
                                DEL_REG+=$CHROM_ARM
                                DEL_REG+=","
                        else
                                ((FOCAL_DEL+=1))
                        fi
		fi
	done < "${TARGET}"

	if [[ $AMP_REG == "" ]]; then
		AMP_REG="NA"
	fi
	if [[ $DEL_REG == "" ]]; then
		DEL_REG="NA"
	fi

        echo -e "$PID\t$LARGE_AMP\t$FOCAL_AMP\t$LARGE_DEL\t$FOCAL_DEL\t$AMP_REG\t$DEL_REG" >> loh.tsv

done < "${FILE}"
