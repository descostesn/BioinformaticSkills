#!/usr/bin/env perl
# This script performs peak calling with hiddendomains software.
# It uses the input file files_wig_X.conf of format: CHROMINFOFILE;BINSIZE;EXPBAM;INPUTBAM;OUTPUTFILEPREFIX
#Nicolas Descostes nov 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
    die("Missing arguments for hiddenDomains_peakCalling.pl\n\n File should contain: CHROMINFOFILE;BINSIZE;EXPBAM;INPUTBAM;OUTPUTFILEPREFIX\n\n");
}

my $chromInfo_file = $arguments_tab[0];
my $bin_size = $arguments_tab[1];
my $expBam = $arguments_tab[2];
my $inputbam = $arguments_tab[3];
my $outputFilePrefix = $arguments_tab[4];

print "This is job number $file_lineNumber \n";

print "/ifs/home/descon01/programs/hiddenDomains/hiddenDomains -g $chromInfo_file -b $bin_size -t $expBam -c $inputbam -o $outputFilePrefix\n\n";

my $commandToLoad = "/ifs/home/descon01/programs/hiddenDomains/hiddenDomains -g $chromInfo_file -b $bin_size -t $expBam -c $inputbam -o $outputFilePrefix";
system($commandToLoad);

