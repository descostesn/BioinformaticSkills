#!/bin/sh

module load meme/4.11.2

perl ./meme-chip.pl $SGE_TASK_ID $NSLOTS $1

echo = `date` job $JOB_NAME done