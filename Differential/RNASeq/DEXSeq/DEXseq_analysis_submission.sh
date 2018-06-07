#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l mem_token=3G
#$ -N DEXSeq
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes jan 2017

# !!-----!! use the following submission command: qsub -pe threaded 1-20 DEXSeq_analysis_submission.sh INPUTFILE!!-----!!   

#This script performs differential alternative splicing by using DEXseq. The Bioconductor package DEXseq implements a method to test for differential exon usage in comparative RNA-Seq experiments.

module load r/3.2.0


perl ./DEXSeq_analysis.pl $SGE_TASK_ID $NSLOTS $1
