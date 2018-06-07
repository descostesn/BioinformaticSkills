#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l mem_token=20G
#$ -N Enrich
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes dec 2016

# !!-----!! use the following submission command: qsub -pe threaded 1-20 chipenrich_analysis_submission.sh INPUTFILE!!-----!!   


# This script performs gene set enrichment analysis from chip-seq peaks using the package Chip-enrich.

module load r/3.2.0

perl ./chipenrich_analysis.pl $SGE_TASK_ID $NSLOTS $1
