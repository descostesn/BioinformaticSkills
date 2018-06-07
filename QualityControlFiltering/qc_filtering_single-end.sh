#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N qc_filtering_single
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# !!-----!! use the following submission command: qsub -pe threaded 1-12 qc_filtering_single-end.sh FILE!!-----!!   

# !!-----!! change path to script !!-----!!

#This script performs the filtering of fastq files according to a defined percentage (X% of sequence) that should have a quality above Y.
# it is calling the perl script qc_filtering_single-end.pl that used the input file list_file_single_qcFiltering.conf


perl ./qc_filtering_single-end.pl $SGE_TASK_ID $NSLOTS $1