#!/usr/bin/env perl
# This script build a markov model using a background file
# It uses the input file files.conf of format: see script]
#Nicolas Descostes March 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 2)
{
    die("Missing arguments for get_background.pl\n\n File should contain: see script\n\n");
}

my $bg_file = $arguments_tab[0];
my $output_file = $arguments_tab[1];

print "This is job number $file_lineNumber \n";

print "fasta-get-markov -m 2 $bg_file > $output_file\n\n";

my $commandToLoad = "fasta-get-markov -m 2 $bg_file > $output_file";
system($commandToLoad);

