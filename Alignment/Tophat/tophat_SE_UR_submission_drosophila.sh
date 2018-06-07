#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N topSEUR
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub -pe threaded 10-20 tophat_SE_UR_submission.sh MISMATCH INPUTFILE !!-----!!   

# !!-----!! change path to script !!-----!!


#This script performs the alignment of single end RNA-seq data, only keeping the uniquely aligned reads
# The argument species will be indicated in the output file name

module load python/2.7
module load bowtie/1.0.0
module load tophat/2.0.9

perl ./tophat_singleEnd_unireads_step1.pl $SGE_TASK_ID $1 $IGENOMES_ROOT/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome $2 $NSLOTS

