#!/usr/bin/env perl
# This script sort bam files
# It uses the input file files_bam_to_sort.conf of format: file_path
#Nicolas Descostes april-june 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


if(scalar(@ARGV) != 3)
{
	die("Missing arguments for sort_bam.pl\n\n\t Command: sort_bam_submission.sh INPUTFILE");	
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

print "sort_bam $file_path $nb_cpu\n";


my $commandToLoad = "/ifs/home/descon01/cluster/scripts/sort_bam_byqueryname $file_path $nb_cpu";
system($commandToLoad);
