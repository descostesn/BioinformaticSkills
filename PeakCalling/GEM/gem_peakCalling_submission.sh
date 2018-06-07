#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l mem_token=30G
#$ -N GEM
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub -pe threaded 1-20 gem_peakCalling_submission.sh INPUTFILE!!-----!!   


#  This script performs peak calling with gem.

module load r/3.3.2
module load java/1.7

perl ./gem_peakCalling.pl $SGE_TASK_ID $NSLOTS $1
