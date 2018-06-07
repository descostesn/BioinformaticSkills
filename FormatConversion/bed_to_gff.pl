#!/usr/bin/env perl
#This script converts a bed file to a gff file
# It uses the input file files_bed_X.conf of format: PATHINPUT;INPUTFILE;THIRDCOLSTRING
#Nicolas Descostes july 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 3)
{
	die("Missing arguments for bed_to_gff.pl\n\n File should contain: PATHINPUT;INPUTFILE;THIRDCOLSTRING\n\n");
}


my $path_to_input = $arguments_tab[0];
my $input_file = $arguments_tab[1];
my $third_col_string = $arguments_tab[2];


print "This is job number $file_lineNumber \n";

print "perl /ifs/home/descon01/cluster/scripts/bed_to_gff_corescript.pl $path_to_input $input_file $third_col_string\n\n";

my $commandToLoad = "perl /ifs/home/descon01/cluster/scripts/bed_to_gff_corescript.pl $path_to_input $input_file $third_col_string";
system($commandToLoad);












