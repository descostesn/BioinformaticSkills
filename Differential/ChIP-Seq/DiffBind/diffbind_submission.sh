#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N Diffbind
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#$ -pe threaded 1-20
#$ -q all.q
#Nicolas Descostes May 2018

# !!-----!! use the following submission command: qsub diffbind_submission.sh INPUTFILE!!-----!!   

#This script performs a differential binding analysis

module load r/3.3.2


perl ./diffbind.pl $SGE_TASK_ID $1
