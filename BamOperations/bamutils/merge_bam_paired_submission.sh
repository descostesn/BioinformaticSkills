#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=60G
#$ -l h_vmem=60G
#$ -l mem_token=60G
#$ -N mergePairedBAM
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub merge_bam_paired_submission.sh INPUTFILE!!-----!!   

#This script merges two bam files representing each pair that were aligned separately

module load ngsutils/0.5.7


perl ./merge_bam_paired.pl $SGE_TASK_ID $1