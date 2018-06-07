#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N Gage_
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes nov 2016

# !!-----!! use the following submission command: qsub gage_analysis_submission.sh INPUTFILE!!-----!!   

#This script performs GSEA analysis with Gage

module load r/3.2.0

perl ./gage_analysis.pl $SGE_TASK_ID $1
