#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N Plot
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes feb 2017

# !!-----!! use the following submission command: qsub seqplots_profiling_submission.sh INPUTFILE!!-----!!   


# This script performs profiling of chip-seq data using annotations.

module load r/3.2.0

perl ./seqplots_profiling.pl $SGE_TASK_ID $1

