#!/usr/bin/env perl
#This script performs the alignment of single end data, only keeping the uniquely aligned reads
# It uses the input file bowtie_single.conf of format: file_path
#Nicolas Descostes april-june 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];
my $species_name = $ARGV[3];
my $mismatch_number = $ARGV[4];
my $bowtie_index_path = $ARGV[5];
my $report_multiple_reads = $ARGV[6];
my $report_unaligned = $ARGV[7];

if(scalar(@ARGV) != 8)
{
	die("Missing arguments for bowtie_singleEnd_unireads.pl\n\n\t Command: bowtie_SE_UR_submission.sh INPUTFILE SPECIES MISMATCH BOWTIEINDEXPATH reportMR(boolean) reportUN(boolean)");	
}

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

#retrive the file name
my @retrieveName = split('/', $file_path);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);
						
print "This is job number $file_lineNumber \n";

print "/ifs/home/descon01/cluster/scripts/UR_bowtie_align $file_path $species_name $nb_cpu $mismatch_number $report_multiple_reads $report_unaligned $bowtie_index_path\n";


my $commandToLoad = "/ifs/home/descon01/cluster/scripts/bash_scripts/UR_bowtie_align $file_path $species_name $nb_cpu $mismatch_number $report_multiple_reads $report_unaligned $bowtie_index_path";
system($commandToLoad);