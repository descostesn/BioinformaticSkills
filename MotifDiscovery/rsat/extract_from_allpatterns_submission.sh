#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l mem_token=3G
#$ -N extract_patterns
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes December 2015

# !!-----!! use the following submission command: qsub extract_from_allpatterns_submission.sh INPUTFILE!!-----!!   


#This script extracts rsat results for each pattern detected

module load r/3.2.0


perl ./extract_from_allpatterns.pl $SGE_TASK_ID $1
