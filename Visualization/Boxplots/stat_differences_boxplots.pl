#!/usr/bin/env perl
# This script performs different statistics (parametric/non parametric) on a table contaning mean expression values that were used to construct boxplots.
# It uses the input file files_peaks_X.conf of format: FILESVALUE;EXPNAMESVEC;OUTPUTFOLDER;REMOVEOUTLIERS
#Nicolas Descostes May 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 4 && scalar(@arguments_tab) != 3)
{
	die("Missing arguments for stat_differences_boxplots.pl\n\n File should contain: FILESVALUE;OUTPUTFOLDER;REMOVEOUTLIERS with optional EXPNAMESVEC\n\n");
}

my $values_file_vec = $arguments_tab[0];
my $output_folder = $arguments_tab[1];
my $remove_outliers = $arguments_tab[2];

print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 4){
	my $expnames_vec = $arguments_tab[4];
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/stat_differences_boxplots.R --valuesFileVec $values_file_vec --expnamesVec $expnames_vec --outputFolder $output_folder --removeOutliers $remove_outliers\n\n";
    my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/stat_differences_boxplots.R --valuesFileVec $values_file_vec --expnamesVec $expnames_vec --outputFolder $output_folder --removeOutliers $remove_outliers";
    system($commandToLoad);
	
}else{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/stat_differences_boxplots.R --valuesFileVec $values_file_vec --outputFolder $output_folder --removeOutliers $remove_outliers\n\n";
    my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/stat_differences_boxplots.R --valuesFileVec $values_file_vec --outputFolder $output_folder --removeOutliers $remove_outliers";
    system($commandToLoad);
}




