#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l mem_token=20G
#$ -N EdgeR
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes october 2016

# !!-----!! use the following submission command: qsub differential_analysis_edgeR_submission.sh INPUTFILE!!-----!!   


# This script performs differential analysis using edgeR package.

module load r/3.2.0


perl ./differential_analysis_edgeR.pl $SGE_TASK_ID $1
