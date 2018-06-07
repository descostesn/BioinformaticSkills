#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N indexBAM
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub create_index_submission.sh INPUTFILE!!-----!!   

#This script creates an index with samtools

module load samtools/1.2.1

perl ./create_index.pl $SGE_TASK_ID $1