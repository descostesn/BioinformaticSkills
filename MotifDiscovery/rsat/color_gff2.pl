#!/usr/bin/env perl
#This script colors sequences according to a motif coordinates given in a gff file
# It uses the input file files_X.conf of format: ANNOFILE;COORDFORCOLORING;PROFILELENGTHBEFORE;PROFILELENGTHAFTER;LOCATION;OUTPUTFOLDER;NAMEEXP
#Nicolas Descostes May 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 7)
{
	die("Missing arguments for color_gff.pl\n\n File should contain: ANNOFILE;COORDFORCOLORING;PROFILELENGTHBEFORE;PROFILELENGTHAFTER;LOCATION;OUTPUTFOLDER;NAMEEXP\n\n");
}


my $annoFile = $arguments_tab[0];
my $coordinatesForColoring = $arguments_tab[1];
my $profileLengthBefore = $arguments_tab[2];
my $profileLengthAfter = $arguments_tab[3];
my $location = $arguments_tab[4];
my $outputFolder = $arguments_tab[5];
my $nameExp = $arguments_tab[6];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/color_gff2.R --annoFile $annoFile --coordinatesForColoring $coordinatesForColoring --profileLengthBefore $profileLengthBefore --profileLengthAfter $profileLengthAfter --location $location --outputFolder $outputFolder --nameExp $nameExp\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/color_gff2.R --annoFile $annoFile --coordinatesForColoring $coordinatesForColoring --profileLengthBefore $profileLengthBefore --profileLengthAfter $profileLengthAfter --location $location --outputFolder $outputFolder --nameExp $nameExp";
system($commandToLoad);





 



