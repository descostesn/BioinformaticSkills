#!/usr/bin/env perl
#This script performs the alignment of single end data using tophat, only keeping the uniquely aligned reads
# It uses the input file tophat_single.conf of format: FILEPATH;EXPNAME
#Nicolas Descostes april-june 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $mismatch_number = $ARGV[1];
my $bowtie_index_path = $ARGV[2];
my $input_file = $ARGV[3];
my $nb_cpu = $ARGV[4];


if(scalar(@ARGV) != 5)
{
	die("Missing arguments for bowtie_singleEnd_unireads.pl\n\n\t Command: tophat_SE_UR_submission.sh MISMATCH INPUTFILE");	
}

#Retrieve the parameters
my $arguments_file = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $arguments_file;

my @arguments_tab = split(';', $arguments_file);

if(scalar(@arguments_tab) != 2)
{
	die("Missing arguments for tophat_singleEnd_unireads_step1.pl\n\n\t File should contain: FILEPATH;EXPNAME");
}


my $file_path = $arguments_tab[0];
my $expname = $arguments_tab[1];

#retrive the file name
my @retrieveName = split('/', $file_path);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);

						
print "This is job number $file_lineNumber \n";

print "tophat -N $mismatch_number --bowtie1 -o $out_path/tophat_singleEnd_UNI/$expname -p $nb_cpu $bowtie_index_path $file_path\n";

my $commandToLoad1 = "mkdir $out_path/tophat_singleEnd_UNI";
system($commandToLoad1);
my $commandToLoad2 = "mkdir $out_path/tophat_singleEnd_UNI/$expname";
system($commandToLoad2); 

my $commandToLoad3 = "tophat -N $mismatch_number --bowtie1 -o $out_path/tophat_singleEnd_UNI/$expname -p $nb_cpu $bowtie_index_path $file_path";
system($commandToLoad3);