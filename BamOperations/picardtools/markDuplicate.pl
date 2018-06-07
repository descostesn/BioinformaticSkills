#!/usr/bin/env perl
#This script marks duplicates in bam file
# It uses the input file files_bam_X.conf of format: BAMFILE;METRICSFILE
#Nicolas Descostes July 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];
my $path_picard = $ARGV[2];

if(scalar(@ARGV) != 3)
{
	die("Missing arguments for markDuplicate.pl\n\n\t Command: markDuplicate_submission.sh INPUTFILE");	
}

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 2)
{
	die("Missing arguments for markDuplicate.pl\n\n File should contain: BAMFILE;METRICSFILE\n\n");
}


my $input_bam_file = $arguments_tab[0];
my $metrics_file = $arguments_tab[1];

my @arg_split = split("\\.bam", $input_bam_file);
my $output_bam_file = join('',$arg_split[0], "_duplicateMarked.bam");


print "This is job number $file_lineNumber\n";

print "java -jar $path_picard/MarkDuplicates.jar INPUT=$input_bam_file OUTPUT=$output_bam_file METRICS_FILE=$metrics_file\n\n";


my $commandToLoad = "java -jar $path_picard/MarkDuplicates.jar INPUT=$input_bam_file OUTPUT=$output_bam_file METRICS_FILE=$metrics_file ASSUME_SORTED=true";
system($commandToLoad);



