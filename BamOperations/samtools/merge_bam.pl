#!/usr/bin/env perl
#This script merge in bam files
# It uses the input file files_bam_X.conf of format: OUTPUTFILE;INPUT1 INPUT2 INPUT3
#Nicolas Descostes July 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 2)
{
	die("Missing arguments for merge_bam.pl\n\n File should contain: OUTPUTFILE;INPUT1 INPUT2 INPUT3...\n\n");
}


my $output_file = $arguments_tab[0];
my $input_file_list = $arguments_tab[1];



print "This is job number $file_lineNumber\n";

print "samtools merge $output_file $input_file_list\n\n";


my $commandToLoad = "samtools merge $output_file $input_file_list";
system($commandToLoad);

