#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N Gquad
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes May 2016

# !!-----!! use the following submission command: qsub -pe threaded 10-20 g-quadruplex_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script colors sequences according to the presence of g-quads

module load r/3.2.0
module load gcc/5.2.0

perl ./g-quadruplex.pl $SGE_TASK_ID $1 $NSLOTS
