#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -t 1-
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l mem_token=20G
#$ -N ClusPro
#$ -M nicolas.descostes@nyumc.org
#$ -m a
#Nicolas Descostes march 2017

# !!-----!! use the following submission command: qsub clusterProfiler_ontologies_submission.sh INPUTFILE!!-----!!   


# This script performs gene ontologies from gene lists with cluster profiler. It enables also the comparison of ontologies between different samples.

module load r/3.3.2

perl ./clusterProfiler_ontologies.pl $SGE_TASK_ID $1
