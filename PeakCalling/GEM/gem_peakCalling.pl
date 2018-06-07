#!/usr/bin/env perl
# This script performs peak calling with GEM (GPS) software.
# It uses the input file files_wig_X.conf of format: EXPERIMENTBED;INPUTBED;GENOMECHROMSIZEFILE;QVALUE;OUTPUTFOLDER;EXPNAME
#Nicolas Descostes march 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 6)
{
	die("Missing arguments for gem_peakCalling.pl\n\n File should contain: EXPERIMENTBED;INPUTBED;GENOMECHROMSIZEFILE;QVALUE;OUTPUTFOLDER;EXPNAME\n\n");
}

my $experiment_bed = $arguments_tab[0];
my $input_bed = $arguments_tab[1];
my $genome_chromsize_file = $arguments_tab[2];
my $q_value = $arguments_tab[3];
my $output_folder = $arguments_tab[4];
my $expname = $arguments_tab[5];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gem_peakCalling.R --experimentBed $experiment_bed --inputBed $input_bed --genomeChromsizeFile $genome_chromsize_file --qValue $q_value --nbCpu $nb_cpu --outputFolder $output_folder --expname $expname\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gem_peakCalling.R --experimentBed $experiment_bed --inputBed $input_bed --genomeChromsizeFile $genome_chromsize_file --qValue $q_value --nbCpu $nb_cpu --outputFolder $output_folder --expname $expname";
system($commandToLoad);



