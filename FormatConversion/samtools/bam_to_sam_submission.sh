#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N BamSam
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes jan 2017

# !!-----!! use the following submission command: qsub bam_to_sam_submission.sh INPUTFILE!!-----!!   

#This script convert bam files to sam files

module load  samtools/0.1.6


perl ./bam_to_sam.pl $SGE_TASK_ID $1
