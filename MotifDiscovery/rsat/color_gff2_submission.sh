#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N color_gff
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes May 2016

# !!-----!! use the following submission command: qsub color_gff_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script colors sequences according to a motif coordinates given in a gff file

module load r/3.2.0


perl ./color_gff2.pl $SGE_TASK_ID $1
