#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N maSigPro
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes may 2017

# !!-----!! use the following submission command: qsub maSigPro_ant_submission.sh INPUTFILE!!-----!!   


# This script performs time course DE analysis of RNA-seq data with the bioconductor package maSigPro.

module load r/3.2.0


perl ./maSigPro_ant.pl $SGE_TASK_ID $1
