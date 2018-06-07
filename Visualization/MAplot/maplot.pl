#!/usr/bin/env perl
# This script performs ma plot from given objects.
# It uses the input file: OBJECTVEC;EXPNAME;OUTPUTFOLDER;FACTOR;LABEL with optional x ticks and y ticks
#Nicolas Descostes march 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5 && scalar(@arguments_tab) != 7)
{
    die("Missing arguments for maplot.pl\n\n File should contain: OBJECTVEC;EXPNAME;OUTPUTFOLDER;FACTOR;LABEL with optional x ticks and y ticks\n\n");
}

my $objectVec = $arguments_tab[0];
my $expnameVec = $arguments_tab[1];
my $outputFolder = $arguments_tab[2];
my $xAxisMultiplyingFactor = $arguments_tab[3];
my $xLabel = $arguments_tab[4];

print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 5)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maplot.R --objectVec $objectVec --expnameVec $expnameVec --outputFolder $outputFolder --xAxisMultiplyingFactor $xAxisMultiplyingFactor --xLabel $xLabel\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maplot.R --objectVec $objectVec --expnameVec $expnameVec --outputFolder $outputFolder --xAxisMultiplyingFactor $xAxisMultiplyingFactor --xLabel $xLabel";
	system($commandToLoad);
}else{
	
	my $xAxisTicks = $arguments_tab[5];
    my $yAxisTicks = $arguments_tab[6];
    
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maplot.R --objectVec $objectVec --expnameVec $expnameVec --outputFolder $outputFolder --xAxisMultiplyingFactor $xAxisMultiplyingFactor --xLabel $xLabel --xAxisTicks $xAxisTicks --yAxisTicks $yAxisTicks\n\n";
    my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maplot.R --objectVec $objectVec --expnameVec $expnameVec --outputFolder $outputFolder --xAxisMultiplyingFactor $xAxisMultiplyingFactor --xLabel $xLabel --xAxisTicks $xAxisTicks --yAxisTicks $yAxisTicks";
    system($commandToLoad);
    
}
	





