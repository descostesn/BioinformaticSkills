#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -l mem_token=4G
#$ -N BinarizeBed
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Dec 2016

# !!-----!! use the following submission command: qsub step3_binarizeBedFiles_submission.sh INPUTFILE!!-----!!   

# After having computed the binarized bed files, this script learn the HMM and output the different states.

module load r/3.2.0
module load java/1.8
module load chromhmm/1.12


perl ./step3_binarizeBedFiles.pl $SGE_TASK_ID $1