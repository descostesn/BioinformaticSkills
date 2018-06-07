#!/usr/bin/env perl
# This script performs peak calling with sicer  
# It uses the input file filesX.conf of format: INPUTDIR;BEDFILE;CONTROLFILE;OUTPUTDIR;SPECIES;REDUNDANCYTHRESHOLD;WINDOWSIZE;FRAGMENTSIZE;EFFECTIVEGENOMEFRACTION;GAPSIZE;FDR
#Nicolas Descostes Jan 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 11)
{
    die("Missing arguments for peakCalling_SICER.pl\n\n File should contain: INPUTDIR;BEDFILE;CONTROLFILE;OUTPUTDIR;SPECIES;REDUNDANCYTHRESHOLD;WINDOWSIZE;FRAGMENTSIZE;EFFECTIVEGENOMEFRACTION;GAPSIZE;FDR\n\n");
}

my $InputDir = $arguments_tab[0];
my $bedFile = $arguments_tab[1];
my $controlFile = $arguments_tab[2];
my $OutputDir = $arguments_tab[3];
my $Species = $arguments_tab[4];
my $redundancyThreshold = $arguments_tab[5];
my $windowSize = $arguments_tab[6];
my $fragmentSize = $arguments_tab[7];
my $effectiveGenomeFraction = $arguments_tab[8];
my $gapSize = $arguments_tab[9];
my $FDR = $arguments_tab[10];

print "This is job number $file_lineNumber \n";

print "bash /ifs/home/descon01/programs/SICER/SICER.sh $InputDir $bedFile $controlFile $OutputDir $Species $redundancyThreshold $windowSize $fragmentSize $effectiveGenomeFraction $gapSize $FDR\n\n";

my $commandToLoad = "bash /ifs/home/descon01/programs/SICER/SICER.sh $InputDir $bedFile $controlFile $OutputDir $Species $redundancyThreshold $windowSize $fragmentSize $effectiveGenomeFraction $gapSize $FDR";
system($commandToLoad);
