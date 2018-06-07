#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N Pipe
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes August 2017

# !!-----!! use the following submission command: qsub pipeline_fastq_to_bigwigs_submission.sh INPUTFILE!!-----!!   


module load r/3.2.0

perl ./pipeline_fastq_to_bigwigs.pl $SGE_TASK_ID $1
