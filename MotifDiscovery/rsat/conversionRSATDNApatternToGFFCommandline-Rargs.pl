#!/usr/bin/env perl
# This script converts rsat results to  gff files
# It uses the input file files_X.conf of format: FOLDERVEC;ANNOFILEVECTOR;OUTPUTFOLDERVEC;LOCATION;PROFILELENGTHBEFORE;PROFILELENGTHAFTER
#Nicolas Descostes November 2015

use strict;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);


if(scalar(@arguments_tab) != 6)
{
	die("Missing arguments for conversionRSATDNApatternToGFFCommandline-Rargs.pl\n\n\t File should contain: FOLDERVEC;ANNOFILEVECTOR;OUTPUTFOLDERVEC;LOCATION;PROFILELENGTHBEFORE;PROFILELENGTHAFTER");
}

my $folderVector = $arguments_tab[0];
my $annoFileVector = $arguments_tab[1];
my $outputFolderVector = $arguments_tab[2];
my $location = $arguments_tab[3];
my $profileLengthBeforeVec = $arguments_tab[4];
my $profileLengthAfterVec = $arguments_tab[5];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/conversionRSATDNApatternToGFFCommandline-Rargs.R --folderVector $folderVector --annoFileVector $annoFileVector --outputFolderVector $outputFolderVector --location $location --profileLengthBeforeVec $profileLengthBeforeVec --profileLengthAfterVec $profileLengthAfterVec\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/conversionRSATDNApatternToGFFCommandline-Rargs.R --folderVector $folderVector --annoFileVector $annoFileVector --outputFolderVector $outputFolderVector --location $location --profileLengthBeforeVec $profileLengthBeforeVec --profileLengthAfterVec $profileLengthAfterVec";
system($commandToLoad); 


	
