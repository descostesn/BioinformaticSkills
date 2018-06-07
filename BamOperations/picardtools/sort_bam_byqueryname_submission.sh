#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=35G
#$ -l h_vmem=35G
#$ -l mem_token=35G
#$ -N sortBAM
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub -pe threaded 10-20 sort_bam_submission.sh INPUTFILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script sort bam files

module load picard-tools/1.88


perl /ifs/home/descon01/scripts_to_submit/sort_bam/sort_bam_byqueryname.pl $SGE_TASK_ID $NSLOTS $1
