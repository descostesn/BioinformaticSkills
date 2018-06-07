#!/usr/bin/env perl
#This script performs union on bed files of intervals
# It uses the input file_X.conf of format: BEDFILEVEC;OUTPUTPATH;OUTPUTFILE
#Nicolas Descostes February 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 4)
{
	die("Missing arguments for makeUnion_bedFiles.pl\n\n File should contain: BEDFILEVEC;OUTPUTPATH;OUTPUTFILE;INPUTFORMAT\n\n");
}

my $file_vec = $arguments_tab[0];
my $output_path = $arguments_tab[1];
my $output_file = $arguments_tab[2];
my $input_format = $arguments_tab[3];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/makeUnion.R --FilesPathVec $file_vec --outputPath $output_path --outputFile $output_file --inputFormat $input_format\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/makeUnion.R --FilesPathVec $file_vec --outputPath $output_path --outputFile $output_file --inputFormat $input_format";
system($commandToLoad);


