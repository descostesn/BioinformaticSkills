#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l mem_token=3G
#$ -N STAR-ant
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes November 2015


# !!-----!! use the following submission command: qsub -pe threaded 1-20 alignment_STAR_SE_submission.sh INPUTFILE!!-----!!


 module load star/2.5.0c

perl ./alignment_STAR_SE.pl $SGE_TASK_ID $1 $NSLOTS /ifs/home/descon01/programs/STAR_index/Hsal85