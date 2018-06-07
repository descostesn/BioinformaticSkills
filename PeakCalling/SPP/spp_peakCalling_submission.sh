#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N SPP
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub -pe threaded 1-20 spp_peakCalling_submission.sh INPUTFILE!!-----!!   


#  This script performs peak calling with spp.

module load r/3.3.2


perl ./spp_peakCalling.pl $SGE_TASK_ID $NSLOTS $1
