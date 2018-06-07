#!/usr/bin/env perl
# This script enables to create a matrix with deeptools
# .conf file should contain: BIGWIGVECTOR;PROFILELENGTHBEFORE;PROFILELENGTHAFTER;OUTPUTMATRIX;BEDFILE;REFERENCEPOINT;SORTINGORDER;CHOICE with option OUTPUTBEDFILE
# Descostes April 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 7 && scalar(@arguments_tab) != 8)
{
    die("Missing arguments for computeMatrix_deeptools.pl\n\n File should contain: BIGWIGVECTOR;PROFILELENGTHBEFORE;PROFILELENGTHAFTER;OUTPUTMATRIX;BEDFILE;REFERENCEPOINT;SORTINGORDER with option OUTPUTBEDFILE\n\n");
}

my $bigwigs_vector = $arguments_tab[0];
my $profileLength_before = $arguments_tab[1];
my $profileLength_after = $arguments_tab[2];
my $output_matrix = $arguments_tab[3];
my $bed_file = $arguments_tab[4];
my $reference_point = $arguments_tab[5];
my $sorting_order = $arguments_tab[6];

print "This is job number $file_lineNumber \n";

#Sorting and outputing a bed file

if(scalar(@arguments_tab) == 8){
	
	my $output_bed_file = $arguments_tab[7];
    print "computeMatrix reference-point -S $bigwigs_vector -R $bed_file -a $profileLength_before -b $profileLength_after -out $output_matrix --outFileSortedRegions $output_bed_file  --referencePoint $reference_point  --binSize 1  --sortRegions $sorting_order --sortUsing mean --averageTypeBins mean --numberOfProcessors $nb_cpu";
    my $commandToLoad = "/ifs/home/descon01/.local/bin/computeMatrix reference-point -S $bigwigs_vector -R $bed_file -a $profileLength_before -b $profileLength_after -out $output_matrix --outFileSortedRegions $output_bed_file  --referencePoint $reference_point  --binSize 1  --sortRegions $sorting_order --sortUsing mean --averageTypeBins mean --numberOfProcessors $nb_cpu";
    system($commandToLoad);
}else{
    #using a sorted bed file
    print "computeMatrix reference-point -S $bigwigs_vector -R $bed_file -a $profileLength_before -b $profileLength_after -out $output_matrix --referencePoint $reference_point --binSize 1 --averageTypeBins mean --numberOfProcessors $nb_cpu";
    my $commandToLoad = "/ifs/home/descon01/.local/bin/computeMatrix reference-point -S $bigwigs_vector -R $bed_file -a $profileLength_before -b $profileLength_after -out $output_matrix --referencePoint $reference_point --binSize 1 --averageTypeBins mean --numberOfProcessors $nb_cpu";
    system($commandToLoad);
}



