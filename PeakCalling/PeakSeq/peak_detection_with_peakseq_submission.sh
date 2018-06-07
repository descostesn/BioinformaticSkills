#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N Peakseq
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub peak_detection_with_peakseq_submission.sh INPUTFILE!!-----!!   


# This script first creates a configuration file and then run peakseq peak detection.

module load r/3.2.0

perl ./peak_detection_with_peakseq.pl $SGE_TASK_ID $1
