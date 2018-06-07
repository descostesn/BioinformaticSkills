#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N StatBoxplots
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Feb 2017

# !!-----!! use the following submission command: qsub stat_differences_boxplots_submission.sh INPUTFILE!!-----!!   


# This script performs different statistics (parametric/non parametric) on a table contaning mean expression values that were used to construct boxplots.

module load r/3.2.0


perl ./stat_differences_boxplots.pl $SGE_TASK_ID $1
