#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N diff_analysis_ants
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes February 2016

# !!-----!! use the following submission command: qsub -pe threaded 1-20 differential_analysis_ants_submission.sh INPUTFILE!!-----!!   

#This script performs a differential analysis for non reference genome

module load r/3.2.0


perl ./differential_analysis_ants.pl $SGE_TASK_ID $NSLOTS $1
