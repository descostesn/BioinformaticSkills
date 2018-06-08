#!/bin/bash

# use the following command: bash rsat-dna-pattern_submission.sh file.conf

#for i in `seq 1 3`
for i in `seq 2 4`
do
	perl ./rsat-dna-pattern.pl $i $1
done