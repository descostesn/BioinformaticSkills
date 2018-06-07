#!/usr/bin/env perl
#This script convert bam files to sam files
# It uses the input file files_sam_X.conf of format: file_path
#Nicolas Descostes jan 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


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

print "samtools view -h -o $out_path/$file_name.sam $file_path\n";

my $commandToLoad = "samtools view -h -o $out_path/$file_name.sam $file_path";
system($commandToLoad);

