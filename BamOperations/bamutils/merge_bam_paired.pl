#!/usr/bin/env perl
#This script merges two bam files representing each pair that were aligned separately
# It uses the input file files_bam_X.conf of format: OUTPUTBAM;BAMFILE1;BAMFILE2;DISCARDBAM1;DISCARDBAM2
#Nicolas Descostes July 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
	die("Missing arguments for merge_bam_paired.pl\n\n File should contain: OUTPUTBAM;BAMFILE1;BAMFILE2;DISCARDBAM1;DISCARDBAM2\n\n");
}

my $output_bam_file = $arguments_tab[0];
my $input_bam_file1 = $arguments_tab[1];
my $input_bam_file2 = $arguments_tab[2];
my $reads_discarded_pair1 = $arguments_tab[3];
my $reads_discarded_pair2 = $arguments_tab[4];

print "This is job number $file_lineNumber\n";

print "bamutils pair -size 50-100000000 -fail1 $reads_discarded_pair1 -fail2 $reads_discarded_pair2 $output_bam_file $input_bam_file1 $input_bam_file2\n\n";


my $commandToLoad = "bamutils pair -size 50-100000000 -fail1 $reads_discarded_pair1 -fail2 $reads_discarded_pair2 $output_bam_file $input_bam_file1 $input_bam_file2";
system($commandToLoad);




