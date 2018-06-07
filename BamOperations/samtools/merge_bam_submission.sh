#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=70G
#$ -l h_vmem=70G
#$ -l mem_token=70G
#$ -N mergeBam
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub merge_bam_submission.sh INPUTFILE!!-----!!   

# This script merge bam files with samtools

module load samtools/1.2.1

perl ./merge_bam.pl $SGE_TASK_ID $1