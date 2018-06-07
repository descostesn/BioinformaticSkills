#!/usr/bin/env perl
# This script performs peak calling with spp.
# It uses the input file files_wig_X.conf of format: BAMFILE;INPUTFILE;OUTPUTFOLDERFORPEAKS;OUTPUTFOLDERFORWIGS;EXPNAME;MINDISTBETWEENPEAKS with options WINDOWSIZEFORDETECTION;ZSCORETHRESHOLD;FDRVALUE
#Nicolas Descostes march 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 6 || scalar(@arguments_tab) != 9)
{
	die("Missing arguments for spp_peakCalling.pl\n\n File should contain: BAMFILE;INPUTFILE;OUTPUTFOLDERFORPEAKS;OUTPUTFOLDERFORWIGS;EXPNAME;MINDISTBETWEENPEAKS with options WINDOWSIZEFORDETECTION;ZSCORETHRESHOLD;FDRVALUE\n\n");
}

my $bam_file = $arguments_tab[0];
my $input_file_script = $arguments_tab[1];
my $output_folder_for_peaks = $arguments_tab[2];
my $output_folder_for_wigs = $arguments_tab[3];
my $expname = $arguments_tab[4];
my $min_distance_between_peaks = $arguments_tab[5];


print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 6)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/spp_peakCalling.R --nbCpu $nbCPU --bamFile $bam_file --inputFile $input_file_script --outputFolderForPeaks $output_folder_for_peaks --outputFolderForWigs $output_folder_for_wigs --expname $expname --minDistanceBetweenPeaks $min_distance_between_peaks\n\n";

	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/spp_peakCalling.R --nbCpu $nbCPU --bamFile $bam_file --inputFile $input_file_script --outputFolderForPeaks $output_folder_for_peaks --outputFolderForWigs $output_folder_for_wigs --expname $expname --minDistanceBetweenPeaks $min_distance_between_peaks";
	system($commandToLoad);
}else{
	
	my $window_size_forPeakDetection = $arguments_tab[6];
	my $zscore_threshold = $arguments_tab[7];
	my $fdr_value = $arguments_tab[8];
	
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/spp_peakCalling.R --nbCpu $nbCPU --bamFile $bam_file --inputFile $input_file_script --outputFolderForPeaks $output_folder_for_peaks --outputFolderForWigs $output_folder_for_wigs --expname $expname --minDistanceBetweenPeaks $min_distance_between_peaks --windowSizeForPeakDetection $window_size_forPeakDetection --zscoreThreshold $zscore_threshold --fdrValue $fdr_value\n\n";

	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/spp_peakCalling.R --nbCpu $nbCPU --bamFile $bam_file --inputFile $input_file_script --outputFolderForPeaks $output_folder_for_peaks --outputFolderForWigs $output_folder_for_wigs --expname $expname --minDistanceBetweenPeaks $min_distance_between_peaks --windowSizeForPeakDetection $window_size_forPeakDetection --zscoreThreshold $zscore_threshold --fdrValue $fdr_value";
	system($commandToLoad);
}





