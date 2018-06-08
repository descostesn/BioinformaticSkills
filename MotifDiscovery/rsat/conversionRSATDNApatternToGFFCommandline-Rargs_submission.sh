#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=1G
#$ -l h_vmem=1G
#$ -l mem_token=1G
#$ -N rsat_to_gff
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes November 2015

# !!-----!! use the following submission command: qsub conversionRSATDNApatternToGFFCommandline-Rargs_submission.sh INPUTFILE!!-----!!   


#This script converts rsat results to  gff files

module load r/3.2.0


perl ./conversionRSATDNApatternToGFFCommandline-Rargs.pl $SGE_TASK_ID $1
