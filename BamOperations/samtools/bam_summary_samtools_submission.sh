#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l mem_token=2G
#$ -N bamstats
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub bam_summary_samtools_submission.sh INPUTFILE!!-----!!   

#This script performs statistics of alignment on a bam file

module load samtools/1.2.1

perl ./bam_summary_samtools.pl $SGE_TASK_ID $1