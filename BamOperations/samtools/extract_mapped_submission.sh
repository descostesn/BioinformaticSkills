#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l mem_token=2G
#$ -N mappedExtract
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub 2_extract_mapped_submission.sh INPUTFILE!!-----!!   

#This script extract mapped reads from a bam file

module load samtools/1.2.1

perl ./2_extract_mapped.pl $SGE_TASK_ID $1