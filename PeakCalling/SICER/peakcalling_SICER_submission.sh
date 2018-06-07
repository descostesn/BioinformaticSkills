#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-1
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N SICER
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Jan 2018

# !!-----!! use the following submission command: qsub peakcalling_SICER_submission.sh INPUTFILE !!-----!!   

#

module load python/2.7

perl ./peakcalling_SICER.pl $SGE_TASK_ID $1
