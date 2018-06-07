#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N GimmeMotifs
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes feb 2018

# !!-----!! use the following submission command: qsub -pe openmpi 10-20 meme_submission.sh INPUTFILE!!-----!!   

#This script performs motif discovery with meme


module load meme/4.11.2

perl ./meme.pl $SGE_TASK_ID $NSLOTS $1
