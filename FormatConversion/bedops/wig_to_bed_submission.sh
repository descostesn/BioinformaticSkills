#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N WigBed
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub wig_to_bed_submission.sh INPUTFILE!!-----!!   


#This script convert wig files to bed files

module load bedops/v2.4.2 


perl ./wig_to_bed.pl $SGE_TASK_ID $1
