#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N MACSDIFF
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub macs2diff_submission.sh INPUTFILE!!-----!!   

#This script performs differential chip-seq peak analysis with MACS2 software

module load macs/2.1.0.20160215

perl ./macs2diff.pl $SGE_TASK_ID $1