#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N fastq_conversion_array
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes april-june 2015

# 		This script launch the conversion of sra file to fastq files by using the perl script sra_fastq_paired-end.pl. This last one use sratoolkit

module load sratoolkit/2.3.2-5

# !!-----!! change path to script !!-----!!

perl path/sra_fastq_single-end.pl $SGE_TASK_ID $1