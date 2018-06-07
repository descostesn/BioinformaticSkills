#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N Union_
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes February 2016

# !!-----!! use the following submission command: qsub makeUnion_bedFiles_submission.sh INPUTFILE!!-----!!   

#This script performs union on bed files of intervals

module load r/3.2.0

perl ./makeUnion.pl $SGE_TASK_ID $NSLOTS $1
