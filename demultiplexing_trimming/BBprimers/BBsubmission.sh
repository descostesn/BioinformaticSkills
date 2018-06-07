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
#$ -q all.q
#Nicolas Descostes april-june 2015

# Submission command: qsub demultiplexing_submission_BB.sh INPUTFILE   
#This script demultiplex a fastq file

perl ./demultiplexing_BB.pl $SGE_TASK_ID $1
