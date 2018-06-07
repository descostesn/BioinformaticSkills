#!/usr/bin/env perl
# This script performs gene ontologies from gene lists with cluster profiler. It enables also the comparison of ontologies between different samples.
# It uses the input file files.conf of format: GFFVEC;EXPNAMESVEC;SPECIESNAME;OUTPUTFOLDER;OUTPUTFORMAT
#Nicolas Descostes march 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
	die("Missing arguments for clusterProfiler_ontologies.pl\n\n File should contain: GFFVEC;EXPNAMESVEC;SPECIESNAME;OUTPUTFOLDER;OUTPUTFORMAT\n\n");
}

my $gff_vec = $arguments_tab[0];
my $expnames_vec = $arguments_tab[1];
my $species_name = $arguments_tab[2];
my $output_folder = $arguments_tab[3];
my $output_format = $arguments_tab[4];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/clusterProfiler_ontologies.R --gffVec $gff_vec --expnamesVec $expnames_vec --speciesName $species_name --outputFolder $output_folder --outputFormat $output_format\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/clusterProfiler_ontologies.R --gffVec $gff_vec --expnamesVec $expnames_vec --speciesName $species_name --outputFolder $output_folder --outputFormat $output_format";
system($commandToLoad);






