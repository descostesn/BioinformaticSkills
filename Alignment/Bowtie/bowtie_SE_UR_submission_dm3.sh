#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-14
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N bwSEUR

# !!-----!! use the following submission command: qsub -pe threaded 10-20 bowtie_SE_UR_submission.sh INPUTFILE SPECIES MISMATCH !!-----!!   

# !!-----!! change path to script !!-----!!

# the "false false" correspond to  "reportMR reportUN" and can be changed

#This script performs the alignment of single end data, only keeping the uniquely aligned reads
# The argument species will be indicated in the output file name

module load  bowtie/1.0.0
module load igenomes/1

perl ./bowtie_singleEnd_unireads.pl $SGE_TASK_ID $NSLOTS $1 $2 $3 $IGENOMES_ROOT/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome false false