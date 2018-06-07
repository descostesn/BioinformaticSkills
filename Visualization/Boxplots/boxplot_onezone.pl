#!/usr/bin/env perl
# This script performs violin plots with the package NGSprofiling.
# It uses the input file files_wig_X.conf of format: see list of args
#Nicolas Descostes feb 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];



#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 10)
{
    die("Missing arguments for boxplots_withggplots.pl\n\n see list of args \n\n");
}


my $bigwigsVec = $arguments_tab[0];
my $expnamesVec = $arguments_tab[1];
my $expnamesRatioVec = $arguments_tab[2];
my $colorsTab = $arguments_tab[3];
my $mainTitle = $arguments_tab[4];
my $outputFolder = $arguments_tab[5];
my $yLabelVec = $arguments_tab[6];
my $yLabelRatioVec = $arguments_tab[7];
my $gffVec = $arguments_tab[8];
my $nameOutputFileNoextVec = $arguments_tab[9];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/boxplot_onezone.R --bigwigsVec $bigwigsVec --expnamesVec $expnamesVec --expnamesRatioVec $expnamesRatioVec --colorsTab $colorsTab --mainTitle $mainTitle --outputFolder $outputFolder --yLabelVec $yLabelVec --yLabelRatioVec $yLabelRatioVec --gffVec $gffVec --nameOutputFileNoextVec $nameOutputFileNoextVec\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/boxplot_onezone.R --bigwigsVec $bigwigsVec --expnamesVec $expnamesVec --expnamesRatioVec $expnamesRatioVec --colorsTab $colorsTab --mainTitle $mainTitle --outputFolder $outputFolder --yLabelVec $yLabelVec --yLabelRatioVec $yLabelRatioVec --gffVec $gffVec --nameOutputFileNoextVec $nameOutputFileNoextVec";
system($commandToLoad);





