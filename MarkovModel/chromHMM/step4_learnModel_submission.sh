#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -l mem_token=4G
#$ -N LearnModel
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Dec 2016

# !!-----!! use the following submission command: qsub -pe threaded 1-20 step4_learnModel_submission.sh INPUTFILE!!-----!!   

# 

module load r/3.2.0
module load java/1.8
module load chromhmm/1.12


perl ./step4_learnModel.pl $SGE_TASK_ID $NSLOTS $1