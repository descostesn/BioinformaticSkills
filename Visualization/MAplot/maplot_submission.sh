#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N MAPlot
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2018

# !!-----!! use the following submission command: qsub maplot_submission.sh INPUTFILE!!-----!!   

module load r/3.2.0

perl ./maplot.pl $SGE_TASK_ID $1

