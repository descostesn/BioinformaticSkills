#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N MAtransform
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes August 2016

# !!-----!! use the following submission command: qsub MAnormPython_diffBinding_submission.sh INPUTFILE!!-----!!   


# This script perform differential binding analysis on chip-seq data using MA transform

module load r/3.2.0
module load python/2.7.3

perl ./MAnormPython_diffBinding.pl $SGE_TASK_ID $1
