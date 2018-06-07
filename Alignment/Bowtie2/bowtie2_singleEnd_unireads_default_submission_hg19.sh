#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l mem_token=10G
#$ -N bow2SEUR
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes october 2015

# !!-----!! use the following submission command: qsub -pe threaded 1-20 bowtie2_singleEnd_unireads_default_submission.sh INPUTFILE SPECIES!!-----!!   

#This script performs the alignment of single end data with bowtie2, only using the default parameters
# The argument species will be indicated in the output file name

module load bowtie2/2.1.0
module load igenomes/1

perl ./bowtie2_singleEnd_unireads_default.pl $SGE_TASK_ID $NSLOTS $1 $2 $IGENOMES_ROOT/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome