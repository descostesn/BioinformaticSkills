#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l mem_token=30G
#$ -N bowSEUR
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub -pe threaded 10-20 bowtie_SE_UR_submission.sh INPUTFILE SPECIES MISMATCH!!-----!!   

# !!-----!! change path to script !!-----!!

# the "false false" correspond to  "reportMR reportUN" and can be changed

#This script performs the alignment of single end data, only keeping the uniquely aligned reads
# The argument species will be indicated in the output file name

module load  bowtie/1.0.0

perl ./bowtie_singleEnd_unireads.pl $SGE_TASK_ID $NSLOTS $1 $2 $3 /ifs/home/descon01/programs/bowtie_index_ant/Hsal85/Hsal85 false false
