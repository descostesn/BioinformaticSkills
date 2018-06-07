#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -l mem_token=4G
#$ -N Analysis
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Dec 2016

# !!-----!! use the following submission command: qsub model_analysis_submission.sh INPUTFILE!!-----!!   

# This script performs different analysis and statistics on the defined markov model.

module load r/3.2.0


perl ./model_analysis.pl $SGE_TASK_ID $1
