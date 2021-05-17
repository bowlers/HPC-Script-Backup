#!/bin/bash
#SBATCH --job-name=delete
#SBATCH --partition=shared
##SBATCH --partition=bioinfo
##SBATCH --account=bioinfo
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --error=./logs/del-%A_%a.err
#SBATCH --output=./logs/del-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

x=1
for d in */ ; do
	if [ "${x}" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
		iName="${d%?}"
		break
	fi
	((x++))
done

rm -r "${iName}"
