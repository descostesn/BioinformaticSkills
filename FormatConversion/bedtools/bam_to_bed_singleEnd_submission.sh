#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N Bam2bed
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes June 2016

# !!-----!! use the following submission command: qsub bam_to_bed_singleEnd_submission.sh INPUTFILE!!-----!!   


#This script converts a bam file to a bed file

module load bedtools/2.22.0

perl ./bam_to_bed_singleEnd.pl $SGE_TASK_ID $1
