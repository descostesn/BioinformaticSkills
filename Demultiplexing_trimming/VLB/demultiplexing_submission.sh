#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N demultiplex
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub demultiplexing_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script demultiplex a fastq file


perl path/demultiplexing.pl $SGE_TASK_ID $1
