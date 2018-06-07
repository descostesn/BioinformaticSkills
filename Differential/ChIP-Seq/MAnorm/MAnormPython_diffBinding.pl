#!/usr/bin/env perl
#This script perform differential binding analysis on chip-seq data using MA transform
# It uses the input file files_X.conf of format: NUMERATORPEAKS;DENOMINATORPEAKS;NUMERATORREADS;DENOMINATORREADS;EXTENSIONNUMERATOR;EXTENSIONDENOMINATOR;NBRANDOMMUTATION;OUTPUTNAME;EXTENSIONMEAN;SUMMITSUMMITDISTANCE;PVALUE;MVALUEBIASED;MVALUEUNBIASED
#Nicolas Descostes August 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];
my $nbcpu = $ARGV[2];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 13)
{
	die("Missing arguments for MAnormPython_diffBinding.pl\n\n File should contain: NUMERATORPEAKS;DENOMINATORPEAKS;NUMERATORREADS;DENOMINATORREADS;EXTENSIONNUMERATOR;EXTENSIONDENOMINATOR;NBRANDOMMUTATION;OUTPUTNAME;EXTENSIONMEAN;SUMMITSUMMITDISTANCE;PVALUE;MVALUEBIASED;MVALUEUNBIASED\n\n");
}



my $numerator_peaks = $arguments_tab[0];     				# File path, It should contain at least three columns, which are chromosome name, start and end position of each peak. If the fourth column exists, it should be summit position (relative to peak start). Otherwise, MAnorm will use the center of peaks instead.
my $denominator_peaks = $arguments_tab[1];
my $numerator_reads = $arguments_tab[2];     				# It should be of bed format, in which the first, second, third and sixth columns are the chromosome, start, end and strand, respectively.
my $denominator_reads = $arguments_tab[3];
my $extension_reads_numerator = $arguments_tab[4];  		# which should be set as the average size of DNA fragments, default=100.
my $extension_reads_denominator = $arguments_tab[5];
my $number_random_permutations = $arguments_tab[6];      	# to test the enrichment of overlapping between two peak sets, default=5.
my $output_name = $arguments_tab[7];                    	# Name of this comparison, which will be also used as the name of folder created to store.
my $extension_for_mean = $arguments_tab[8];              	# default=1000, 2*extension=size of the window centered at peak summit to calculate reads density. The window size should match the typical length of peaks, thus we recommend extension=1000 for sharp histone marks like H3K4me2/3 or H3K9/27ac, extension=500 for transcription factor or DNase-seq.
my $summit_to_summit_distance_cutoff = $arguments_tab[9];   # default=extension/2. Only those common peaks with distance between their summits in 2 samples smaller than this value will be considered as real common peaks for building the normalization model.  
my $p_value = $arguments_tab[10];     						# Cutoff of P-value to define biased (high-confidence sample 1 or 2-specific) peaks, default=0.01.
my $m_value_biased = $arguments_tab[11];     				# Cutoff of M-value to define biased peaks, default=1. Sample 1 biased peaks are defined as sample 1 unique peaks with M-value > mcut_biased and P-value < pcut_biased, while sample 2 biased peaks are defined as sample 2 unique peaks with M-value < -1*mcut_biased and P-value < pcut_biased.
my $m_value_unbiased = $arguments_tab[12];   				# Cutoff of M-value to define unbiased (high-confidence non-specific) peaksbetween 2 samples, default=1. They are defined to be the common peaks with -1*mcut_unbiased < M-value < mcut_unbiased and P-value > pcut_biased.



print "This is job number $file_lineNumber \n";

print "python /ifs/home/descon01/cluster/scripts/python_scripts/MAnorm.py --p1 $numerator_peaks --p2 $denominator_peaks --r1 $numerator_reads --r2 $denominator_reads --s1 $extension_reads_numerator --s2 $extension_reads_denominator -n $number_random_permutations -o $output_name -e $extension_for_mean -d $summit_to_summit_distance_cutoff -p $p_value -m $m_value_biased -u $m_value_unbiased\n\n";

my $commandToLoad = "python /ifs/home/descon01/cluster/scripts/python_scripts/MAnorm.py --p1 $numerator_peaks --p2 $denominator_peaks --r1 $numerator_reads --r2 $denominator_reads --s1 $extension_reads_numerator --s2 $extension_reads_denominator -n $number_random_permutations -o $output_name -e $extension_for_mean -d $summit_to_summit_distance_cutoff -p $p_value -m $m_value_biased -u $m_value_unbiased";
system($commandToLoad);




