#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N Boxplots
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes nov 2016

# !!-----!! use the following submission command: qsub -pe threaded 1-20 boxplots_with_rnaseq_submission.sh INPUTFILE!!-----!!   

#This script generates boxplots for rna-seq experiments.

module load r/3.2.0


perl ./boxplots_with_rnaseq.pl $SGE_TASK_ID $NSLOTS $1
