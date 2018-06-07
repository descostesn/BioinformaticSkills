#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N bed2GFF
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes july 2015

# !!-----!! use the following submission command: qsub bed_to_gff_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script converts a bed file to a gff file

module load r/3.2.0

perl ./bed_to_gff.pl $SGE_TASK_ID $1
