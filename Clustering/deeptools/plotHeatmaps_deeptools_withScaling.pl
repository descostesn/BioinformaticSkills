#!/usr/bin/env perl
# This script enables to create heatmaps with the matrix generated with the deeptools script computeMatrix
# .conf file should contain: MATRIXFILE;OUTPUTPDF;SORTREGIONS;COLOR;XAXISLAB;REFPOINTLAB;SAMPLESLABELVEC;PLOTTITLE;ZMIN;ZMAX;YMIN;YMAX
# Descostes April 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 12)
{
    die("Missing arguments for plotHeatmaps_deeptools.pl\n\n File should contain: MATRIXFILE;OUTPUTPDF;SORTREGIONS;COLOR;XAXISLAB;REFPOINTLAB;SAMPLESLABELVEC;PLOTTITLE;ZMIN;ZMAX;YMIN;YMAX\n\n");
}

my $matrix_file = $arguments_tab[0];
my $output_pdf = $arguments_tab[1];
my $sort_regions = $arguments_tab[2];
my $color = $arguments_tab[3];
my $xaxis_lab = $arguments_tab[4];
my $refpoint_label = $arguments_tab[5];
my $samples_label_vector = $arguments_tab[6];
my $plot_title = $arguments_tab[7];
my $ZMIN = $arguments_tab[8];
my $ZMAX = $arguments_tab[9];
my $YMIN = $arguments_tab[10];
my $YMAX = $arguments_tab[11];



print "This is job number $file_lineNumber \n";

print "/ifs/home/descon01/.local/bin/plotHeatmap -m $matrix_file -out $output_pdf --sortRegions $sort_regions --colorMap $color --xAxisLabel $xaxis_lab --refPointLabel $refpoint_label --samplesLabel $samples_label_vector --plotTitle $plot_title --zMin $ZMIN --zMax $ZMAX --yMin $YMIN --yMax $YMAX";
my $commandToLoad = "/ifs/home/descon01/.local/bin/plotHeatmap -m $matrix_file -out $output_pdf --sortRegions $sort_regions --colorMap $color --xAxisLabel $xaxis_lab --refPointLabel $refpoint_label --samplesLabel $samples_label_vector --plotTitle $plot_title --zMin $ZMIN --zMax $ZMAX --yMin $YMIN --yMax $YMAX";
system($commandToLoad);









