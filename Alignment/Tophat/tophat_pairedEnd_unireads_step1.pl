#!/usr/bin/env perl
#This script performs the alignment of paired end data using tophat, only keeping the uniquely aligned reads
# It uses the input file tophat_single.conf of format: file_path1;file_path2;expname
#Nicolas Descostes april-june 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $mismatch_number = $ARGV[1];
my $bowtie_index_path = $ARGV[2];
my $input_file = $ARGV[3];
my $nb_cpu = $ARGV[4];


if(scalar(@ARGV) != 5)
{
	die("Missing arguments for bowtie_pairedEnd_unireads.pl\n\n\t Command: tophat_PE_UR_submission.sh MISMATCH INPUTFILE");	
}


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @paired_files = split(';', $file_path);
my $file_path1 = $paired_files[0];
my $file_path2 = $paired_files[1];
my $expname = $paired_files[2];

#retrive the file name
my @retrieveName = split('/', $file_path1);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);

						
print "This is job number $file_lineNumber \n";

print "tophat -N $mismatch_number --bowtie1 -o $out_path/tophat_pairedEnd/$expname -p $nb_cpu $bowtie_index_path $file_path1,$file_path2\n";


my $commandToLoad = "tophat -N $mismatch_number --bowtie1 -o $out_path/tophat_pairedEnd/$expname -p $nb_cpu $bowtie_index_path $file_path1 $file_path2";
system($commandToLoad);