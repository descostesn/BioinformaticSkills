#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N var2fix
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes July 2015

# !!-----!! use the following submission command: qsub wigVariable2Fixed_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script converts a wig file in variable steps to one in fixed steps

module load r/3.2.0


perl ./wigVariable2Fixed.pl $SGE_TASK_ID $1
