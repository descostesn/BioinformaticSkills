#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N overlapGFF
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes feb 2016

# !!-----!! use the following submission command: qsub vennDiagram_overlapGFF_submission.sh INPUTFILE!!-----!!   

#This script performs overlap between 2 to 4 GFF files

module load r/3.2.0

perl ./vennDiagram_overlapGFF.pl $SGE_TASK_ID $1