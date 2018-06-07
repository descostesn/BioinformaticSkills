#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l mem_token=5G
#$ -N PrePeakseq
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub preprocess_bam_files_PeakSeq_submission.sh INPUTFILE!!-----!!   


# Preprocess bam files into binary files used by PeakSeq by pipping with samtools.

module load samtools/1.3

perl ./preprocess_bam_files_PeakSeq.pl $SGE_TASK_ID $1
