#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=45G
#$ -l h_vmem=45G
#$ -l mem_token=45G
#$ -N Heatmap
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes feb 2017

# !!-----!! use the following submission command: qsub kmean_with_seqplots_submission.sh INPUTFILE!!-----!!   


# This script performs heatmap of kmean clustering of chip-seq data using annotations.

module load r/3.2.0

perl ./kmean_with_seqplots.pl $SGE_TASK_ID $1

