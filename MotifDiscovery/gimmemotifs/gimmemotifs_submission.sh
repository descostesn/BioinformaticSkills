#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N GimmeMotifs
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes dec 2016

# !!-----!! use the following submission command: qsub -pe threaded 15 gimmemotifs_submission.sh INPUTFILE!!-----!!   

#This script performs motif discovery with gimmemotifs


module load metaseq/0.5.5.4
module load meme/4.9.1

perl ./gimmemotifs.pl $SGE_TASK_ID $1
