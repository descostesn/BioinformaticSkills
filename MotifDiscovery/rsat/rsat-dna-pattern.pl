#!/usr/bin/env perl
# This script looks for pattern in fasta sequences file
# It uses the input file files_X.conf of format: FASTAFILE;OUTPUTFILE;PATTERNFILE
#Nicolas Descostes November

use strict;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);


if(scalar(@arguments_tab) != 3)
{
	die("Missing arguments for rsat-dna-pattern.pl\n\n\t File should contain: FASTAFILE;OUTPUTFILE;PATTERNFILE");
}

my $fasta_file = $arguments_tab[0];
my $output_file = $arguments_tab[1];
my $pattern_file = $arguments_tab[2];


print "This is job number $file_lineNumber \n";

print "dna-pattern -i $fasta_file -format fasta -o $output_file -pl $pattern_file -1str -origin 0\n\n";


my $commandToLoad = "dna-pattern -i $fasta_file -format fasta -o $output_file -pl $pattern_file -1str -origin 0";
system($commandToLoad);





