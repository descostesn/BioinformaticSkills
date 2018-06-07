#!/usr/bin/env perl
#This script performs the filtering of fastq files according to a defined percentage (X% of sequence) that should have a quality above Y.
# It uses the input file list_file_single_qcFiltering.conf of format: file_path;percent;threshold
#Nicolas Descostes april-june 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $file_parameters = $ARGV[2];


#Retrieve the parameters
my $parameters = `head -n $file_lineNumber $file_parameters | tail -n1`;
my @parameters = split(";", $parameters);

if(scalar(@parameters) != 3){
	
		die "One or more arguments are missing, end of program \n\n";
}		


my $file_path = $parameters[0];
my $percent = $parameters[1];
my $threshold = $parameters[2];


#retrive the file name
my @retrieveName = split('/', $file_path);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);
						
print "This is job number $file_lineNumber \n";

print " perl /ifs/home/descon01/programs/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -se $file_path N A -l $percent -s $threshold -p $nb_cpu -o $out_path \n";

# Here are the following options:
# -se: N means to not filter for primer/adapter and A stands for automatic detection of fastq variants
# -l:  The cut-off value for percentage of read length that should be of given quality (between 0 and 100)
# -s:  The cut-off value for PHRED quality score for high-quality filtering (between 0 and 40)

my $commandToLoad = "perl /ifs/home/descon01/programs/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -se $file_path N A -l $percent -s $threshold -p $nb_cpu -o $out_path";
system($commandToLoad);