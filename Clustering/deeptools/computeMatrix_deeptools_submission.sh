#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N Deeptools_heatmap
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#$ -q all.q
#$ -pe threaded 10-20
#Nicolas Descostes April 2018

# !!-----!! use the following submission command: qsub computeMatrix_deeptools_submission.sh INPUTFILE!!-----!!   

module load python/3.5.3


perl ./computeMatrix_deeptools.pl $SGE_TASK_ID $NSLOTS $1
