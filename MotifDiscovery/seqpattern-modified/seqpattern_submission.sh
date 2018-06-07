#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l mem_token=20G
#$ -N Seqpattern
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes january 2017

# !!-----!! use the following submission command: qsub -pe threaded 1-20 seqpattern_submission.sh INPUTFILE!!-----!!   


#This script generate heatmaps and av profiles of motifs from a fasta file.

module load r/3.2.0


perl ./seqpattern.pl $SGE_TASK_ID $NSLOTS $1
