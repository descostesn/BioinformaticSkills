#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=100G
#$ -l h_vmem=100G
#$ -l mem_token=100G
#$ -N duplicate
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub markDuplicate_submission.sh INPUTFILE!!-----!!   

# This script mark duplicated reads in a bam file with picard tool

module load picard-tools/1.88

perl ./markDuplicate.pl $SGE_TASK_ID $1 $PICARD_DIR