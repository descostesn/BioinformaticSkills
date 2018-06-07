#!/usr/bin/env perl
#This script extract mapped reads from a bam file
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
	die("Missing arguments for 2_extract_mapped.pl\n\n File should contain: BAMFILE\n\n");
}


my $input_bam_file = $arguments_tab[0];

my @arg_split = split("\\.bam", $input_bam_file);
my $output_bam_file = join('',$arg_split[0], "_mapped.bam");


print "This is job number $file_lineNumber\n";

print "samtools view -b -F 4 $input_bam_file > $output_bam_file\n\n";


my $commandToLoad = "samtools view -b -F 4 $input_bam_file > $output_bam_file";
system($commandToLoad);
