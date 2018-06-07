#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l mem_token=30G
#$ -N Expression
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes feb 2017

# !!-----!! use the following submission command: qsub multiple_techniques_submission.sh INPUTFILE!!-----!!   


# This script plots the distribution of gene expression from a bam file. define different classes with different classification techniques, and output bed files for each category.

module load r/3.2.0

perl ./multiple_techniques.pl $SGE_TASK_ID $1

