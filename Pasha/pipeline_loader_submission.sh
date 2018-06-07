#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=130G
#$ -l h_vmem=130G
#$ -l mem_token=130G
#$ -N pasha
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub pipeline_loader_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script submits bam files to the pipeline

module load r/3.2.0


perl ./pipeline_loader.pl $SGE_TASK_ID $1
