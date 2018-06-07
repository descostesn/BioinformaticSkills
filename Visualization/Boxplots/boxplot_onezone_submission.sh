#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N Boxplot
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#$ -q gpu0.q
#Nicolas Descostes feb 2018

# !!-----!! use the following submission command: qsub -pe boxplots_onezone_submission.sh INPUTFILE!!-----!!   


module load r/3.4.3-cairo


perl ./boxplot_onezone.pl $SGE_TASK_ID $1
