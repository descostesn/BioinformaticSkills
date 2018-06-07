#!/usr/bin/env perl
#This script converts a bam file to a bed file
# It uses the input file files_X.conf of format: BAMFILE;OUTPUTBED
#Nicolas Descostes june 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 2)
{
	die("Missing arguments for bam_to_bed_singleEnd.pl\n\n File should contain: BAMFILE;OUTPUTBED\n\n");
}


my $bam_files = $arguments_tab[0];
my $output_bed = $arguments_tab[1];

print "This is job number $file_lineNumber \n";

print "bedtools bamtobed -i $bam_files > $output_bed\n";

my $commandToLoad = "bedtools bamtobed -i $bam_files > $output_bed";
system($commandToLoad);
