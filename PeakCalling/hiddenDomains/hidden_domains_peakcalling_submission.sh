#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l mem_token=15G
#$ -N HiddenDomains
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes Nov 2017

# !!-----!! use the following submission command: qsub hidden_domains_peakcalling_submission.sh INPUTFILE!!-----!!   

#This script performs peak detection with hiddendomains software

module load r/3.2.0
module load samtools/1.3
module load perl/5.22.1

perl ./hiddenDomains_peakcalling.pl $SGE_TASK_ID $1