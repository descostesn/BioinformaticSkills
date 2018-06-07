#!/usr/bin/env perl
#This script performs statistics of alignment on a bam file
# It uses the input file files_bam_X.conf of format: BAMFILE
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
	die("Missing arguments for 3_bam_summary_samtools.pl\n\n File should contain: BAMFILE\n\n");
}

my $input_bam_file = $arguments_tab[0];

print "This is job number $file_lineNumber\n";

print "samtools flagstat $input_bam_file\n\n";


my $commandToLoad = "samtools flagstat $input_bam_file";
system($commandToLoad);





