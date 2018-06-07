#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N adapters
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub cutadapt_colorspace_submission.sh INPUTFILE!!-----!!   

#This script removes adapters from color space sequences

module load python/3.4.3

perl ./cutadapt_colorspace.pl $SGE_TASK_ID $1