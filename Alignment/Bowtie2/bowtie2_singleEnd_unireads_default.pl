#!/usr/bin/env perl
#This script performs the alignment of single end data using bowtie2, only using the default parameters
# It uses the input file bowtie2_single.conf of format: FILEPATH
#Nicolas Descostes october 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $nbcpu = $ARGV[1];
my $input_file = $ARGV[2];
my $species = $ARGV[3];
my $bowtie_index = $ARGV[4];



if(scalar(@ARGV) != 5)
{
	die("Missing arguments for bowtie_singleEnd_unireads.pl\n\n\t Command: bowtie2_singleEnd_unireads_default_submission.sh INPUTFILE SPECIES");	
}

#Retrieve the parameters
my $arguments_file = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $arguments_file;

my @arguments_tab = split(';', $arguments_file);

if(scalar(@arguments_tab) != 1)
{
	die("Missing arguments for tophat_singleEnd_unireads_step1.pl\n\n\t File should contain: FILEPATH");
}

my $fastq_file = $arguments_tab[0];

#retrive the file name
my @retrieveName = split('/', $fastq_file);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);
my $output_file = "$out_path/$file_name-bw2align-$species-defaultparam.sam";
my $log_file = "$out_path/$file_name-bw2align-$species-defaultparam_LOG.txt";
		
print "bowtie2 --sensitive -p $nbcpu -x $bowtie_index -U $fastq_file -S $output_file 2> $log_file\n";


my $commandToLoad = "bowtie2 --sensitive -p $nbcpu -x $bowtie_index -U $fastq_file -S $output_file 2> $log_file";
system($commandToLoad);
