#!/usr/bin/env perl
#  This script first creates a configuration file and then run peakseq peak detection.
# It uses the input file files_wig_X.conf of format:  
#Nicolas Descostes march 2017

use strict;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 11)
{
	die("Missing arguments for peak_detection_with_peakseq.pl.pl\n\n");
}

my $exp_id = $arguments_tab[0];
my $tag_length = $arguments_tab[1];
my $target_fdr = $arguments_tab[2];
my $number_simulation = $arguments_tab[3];
my $min_interpeak_distance = $arguments_tab[4];
my $mappability_map_file = $arguments_tab[5];
my $chipseq_reads_folder = $arguments_tab[6];
my $input_reads_folder = $arguments_tab[7];
my $q_value = $arguments_tab[8];
my $background_model_mode = $arguments_tab[9];
my $output_folder = $arguments_tab[10];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/peak_detection_with_peakseq.R --expId $exp_id --tagLength $tag_length --targetFdr $target_fdr --numberSimulation $number_simulation --minInterpeakDistance $min_interpeak_distance --mappabilityMapFile $mappability_map_file --chipseqReadsFolder $chipseq_reads_folder --inputReadsFolder $input_reads_folder --Qvalue $q_value --backgroundModelMode $background_model_mode --outputFolder $output_folder\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/peak_detection_with_peakseq.R --expId $exp_id --tagLength $tag_length --targetFdr $target_fdr --numberSimulation $number_simulation --minInterpeakDistance $min_interpeak_distance --mappabilityMapFile $mappability_map_file --chipseqReadsFolder $chipseq_reads_folder --inputReadsFolder $input_reads_folder --Qvalue $q_value --backgroundModelMode $background_model_mode --outputFolder $output_folder";
system($commandToLoad);




