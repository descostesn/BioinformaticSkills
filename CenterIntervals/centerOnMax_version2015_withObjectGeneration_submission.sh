#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N Center
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes March 2018

# !!-----!! use the following submission command: qsub -pe threaded 5-20 centerOnMax_version2015_withObjectGeneration_submission.sh INPUTFILE!!-----!!   


module load r/3.2.0

perl ./centerOnMax_version2015_withObjectGeneration.pl $SGE_TASK_ID $NSLOTS $1
