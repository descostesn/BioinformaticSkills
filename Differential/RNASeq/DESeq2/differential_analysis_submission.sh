#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N diff_analysis
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes February 2016

# !!-----!! use the following submission command: qsub -pe threaded 1-20 differential_analysis_submission.sh INPUTFILE!!-----!!   

#This script performs a differential analysis

module load r/3.2.0


perl ./differential_analysis.pl $SGE_TASK_ID $NSLOTS $1
