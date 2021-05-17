#!/bin/bash
#SBATCH --job-name=parseVCF
##SBATCH --partition=shared
#SBATCH --partition=bioinfo
#SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --mem=58G
#SBATCH --cpus-per-task=1
#SBATCH --error=./logs/parse/parse-%A-%a.err
#SBATCH --output=./logs/parse/parse-%A-%a.out

module load lang/Anaconda3
source activate $HOME/biocore_lts/scott/envir/maf
module load lang/Java
module load lang/Python

if [ -z "${1}" ]; then
	echo "Provide a study"
	exit 1
else
	STUDY="${1}"
fi

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
	SLURM_ARRAY_TASK_ID=1
fi

if [ -z "${2}" ]; then
	ACTION="PARSE"
else
	ACTION="MERGE"
fi

if [ ! -d $HOME/biocore_lts/scott/"${STUDY}"/somatic ]; then
	mkdir $HOME/biocore_lts/scott/"${STUDY}"/somatic
fi

cd $HOME/biocore_lts/scott/"${STUDY}"/somatic

FILE=design.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID} p" ../"${FILE}")
TUMOR=$(echo $LINE | awk '{print $1}')
NORMAL=$(echo $LINE | awk '{print $2}')
PID=$(echo $LINE | awk '{print $3}')

if [ ! -d ./results ]; then
	mkdir results
fi

if [ "${ACTION}" == "PARSE" ]; then
	if [ ! -f ./"${PID}".filtered.snp.Somatic.hc.vcf ]; then
		source activate $HOME/biocore_lts/scott/envir/env
		varscan processSomatic ../varscan/filtered/"${PID}".filtered.snp.vcf
		cp ../varscan/filtered/"${PID}".filtered.snp.Somatic.hc.vcf .
		cp ../varscan/filtered/"${PID}".filtered.snp.Somatic.vcf .
		source activate $HOME/biocore_lts/scott/envir/maf
	fi
	if [ ! -f ./"${PID}".vcf ]; then
		{ cat "${PID}".filtered.snp.Somatic.hc.vcf; sed '1,18d' "${PID}".filtered.snp.Somatic.vcf; } > "${PID}".vcf
	fi
	if [ ! -f ./"${PID}".annot.vcf ]; then
		java -Xmx12g -jar $HOME/biocore_lts/scott/snpEff/snpEff.jar -noStats hg38 "${PID}".vcf > "${PID}".annot.vcf
	fi
	if [ ! -f ./results/"${PID}".table ]; then
		gatk VariantsToTable -V "${PID}".annot.vcf -F CHROM -F POS -F REF -F ALT -F ANN -GF FREQ -GF DP -GF AD -O ./results/"${PID}".table
	fi
	if [ ! -f ./results/"${PID}".mod.table ]; then
		HEADER="GENE\t"
		HEADER+="${PID}\tDP\tAD"
		echo -e "${HEADER}" > ./results/"${PID}".mod.table
		while read AROW; do
			STRING=${AROW:0:5}
			if [ "${STRING}" == "CHROM" ]; then
				continue
			fi
			ANN=$( echo $AROW | awk '{print $5}' )
			IFS='|' read -r -a ARRAY <<< "${ANN}"
			if [ "${ARRAY[1]}" == "synonymous_variant" ]; then
				continue
			fi
			GENE=${ARRAY[3]}
			FREQ=$( echo $AROW | awk '{print $9}' )
			FREQ=${FREQ::-1}
			RD=$( echo $AROW | awk '{print $10}' )
			AD=$( echo $AROW | awk '{print $11}' )
			echo -e "$GENE\t$FREQ\t$RD\t$AD" >> ./results/"${PID}".mod.table
		done < ./results/"${PID}".table
	fi
	if [ ! -f ./results/"${PID}".sorted.table ]; then
		awk -F'\t' '$1!=""' ./results/"${PID}".mod.table > ./results/"${PID}".mod.table.tmp
		mv ./results/"${PID}".mod.table.tmp ./results/"${PID}".mod.table
		(head -n 1 ./results/"${PID}".mod.table && tail -n +3 ./results/"${PID}".mod.table | sort ) > ./results/"${PID}".sorted.table
	fi
else
	if [ ! -f ./results/header.lst ]; then
		cat ./results/*.sorted.table > ./results/total.table
		echo "GENE" > ./results/header.lst
		cut -f1 ./results/total.table | sort | uniq >> ./results/header.lst
		rm ./results/total.table
	fi

	if [ -z "${3}" ]; then
		python $HOME/biocore_lts/scott/scripts/merge.py
	else
		python $HOME/biocore_lts/scott/scripts/merge_probability.py
	fi
fi
