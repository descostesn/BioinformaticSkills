#!/usr/bin/env perl
#This script demultiplex files by using a barcode
# It uses the input file files_multiplex_X.conf of format: file_fastq
#Nicolas Descostes april-june 2015

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

print "perl /ifs/home/descon01/cluster/scripts/vblsplit.sh $file_path\n\n";


my $commandToLoad = "perl /ifs/home/descon01/cluster/scripts/vblsplit.pl $file_path";
system($commandToLoad);