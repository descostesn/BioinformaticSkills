#!/usr/bin/env perl
#This script performs a differential analysis for non-reference genome
# It uses the input file_X.conf of format: BAMFILE;SINGLE;SAMPLEINFOFILE;NAMEOFREF;OUTPUTFOLDER;COMPARISONNAME;REFSEQANNO;EXPNAME_VEC
#Nicolas Descostes February 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 8)
{
	die("Missing arguments for differential_analysis_ants.pl\n\n File should contain: BAMFILE;SINGLE;SAMPLEINFOFILE;NAMEOFREF;OUTPUTFOLDER;COMPARISONNAME;REFSEQANNO;EXPNAME_VEC\n\n");
}

my $bam_files_vec = $arguments_tab[0];
my $single = $arguments_tab[1];
my $samples_info_file = $arguments_tab[2];
my $name_of_reference = $arguments_tab[3];
my $output_folder = $arguments_tab[4];
my $comparison_name = $arguments_tab[5];
my $refseq_anno = $arguments_tab[6];
my $expnames_vec = $arguments_tab[7];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis_ants.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis_ants.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec";
system($commandToLoad);





