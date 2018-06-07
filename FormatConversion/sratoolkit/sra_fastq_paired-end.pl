#!/usr/bin/env perl
#This script enable to convert a sra to fastq taking into account paired-end tags. it uses a file with input data called list_file_paired_end.conf.
# list_file_paired_end.conf format is: sra-file
#Nicolas Descostes april-june 2015

use strict;
my $lineNumber = $ARGV[0];
my $file_parameters = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $lineNumber $file_parameters | tail -n1`;
chomp $file_path;

#retrive the file name
my @retrieveName = split('/', $file_path);
my @file_name_array = split('\\.',$retrieveName[-1]);
my $file_name = $file_name_array[0];

#retrieve the path
my $out_path =  join('/', @retrieveName[0..$#retrieveName-1]);

print "This is job number $lineNumber \n";
print "About to perform $file_path\n";
print "About to write $file_name.fastq\n";
print "To the output folder $out_path\n";

my $cmd = "fastq-dump --split-files --outdir $out_path $file_path";
system($cmd);