#!/usr/bin/env perl
#This script demultiplex files by using a Illumina barcode
# It uses the input file files_multiplex_X.conf of format: file_fastq
#Nicolas Descostes jan 2017

use strict;
 
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

if(scalar(@ARGV) != 2)
{
    die("Missing arguments for demultiplexing.pl\n\n\t Command: demultiplexing_submission.sh INPUTFILE");   
}

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

print "This is job number $file_lineNumber \n";

print "perl IlluminaSplit.pl $file_path\n\n";


my $commandToLoad = "perl IlluminaSplit.pl $file_path";
system($commandToLoad);