#!/usr/bin/env perl
#This script creates an index with samtools
# It uses the input file files_bam_X.conf of format: BAMFILEPATH
#Nicolas Descostes July 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 1)
{
	die("Missing arguments for create_index.pl\n\n File should contain: BAMFILEPATH\n\n");
}

my $input_bam = $arguments_tab[0];

print "This is job number $file_lineNumber\n";

print "samtools index $input_bam\n\n";
my $commandToLoad2 = "samtools index $input_bam";
system($commandToLoad2);
