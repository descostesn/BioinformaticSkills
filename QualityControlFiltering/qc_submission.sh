#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l mem_token=2G
#$ -N QC_array
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# command: qsub qc_submission.sh FILE

#This script is loading the generation of a fastqc report in order to filter the reads according to their quality of sequencing. It uses the fastqc toolkit and is calling the quality_control.pl script.

module load fastqc/0.11.4

# !!-----!! change path to script !!-----!!

perl path/quality_control.pl $SGE_TASK_ID $1