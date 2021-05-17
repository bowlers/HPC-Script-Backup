#!/bin/bash
#SBATCH --job-name=bt2build
#SBATCH --partition=shared
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --error=./logs/del-%A_%a.err
#SBATCH --output=./logs/del-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sbowler@hawaii.edu

module load bio/Bowtie2

cd $HOME/biocore_lts/scott/resources

bowtie2-build GRCh38.fa GRCh38
