#!/usr/bin/env perl
#  This script performs differential alternative splicing by using DEXseq. The Bioconductor package DEXseq implements a method to test for differential exon usage in comparative RNA-Seq experiments.
# It uses the input file files_correl_X.conf of format: COUNTFILEVEC;FLATTENEDFILE;ROWNAMEVEC;CONDITIONVEC;LIBRARYTYPEVEC;OUTPUTFOLDER;FDRCUTOFF 
#Nicolas Descostes jan 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $nbcpu  = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 7)
{
	die("Missing arguments for DEXSeq_analysis.pl\n\n File should contain: COUNTFILEVEC;FLATTENEDFILE;ROWNAMEVEC;CONDITIONVEC;LIBRARYTYPEVEC;OUTPUTFOLDER;FDRCUTOFF\n\n");
}


my $count_file_vec = $arguments_tab[0];
my $flattened_file = $arguments_tab[1];
my $row_name_vec = $arguments_tab[2];
my $condition_vec = $arguments_tab[3];
my $library_type_vec = $arguments_tab[4];
my $output_folder = $arguments_tab[5];
my $fdr_cutoff = $arguments_tab[6];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/DEXseq_analysis.R --countFilesVec $count_file_vec --flattenedFile $flattened_file --rowNamesVec $row_name_vec --conditionVec $condition_vec --libtypeVec $library_type_vec --outputFolder $output_folder --fdrCutoff $fdr_cutoff --nbCpu $nbcpu\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/DEXseq_analysis.R --countFilesVec $count_file_vec --flattenedFile $flattened_file --rowNamesVec $row_name_vec --conditionVec $condition_vec --libtypeVec $library_type_vec --outputFolder $output_folder --fdrCutoff $fdr_cutoff --nbCpu $nbcpu";
system($commandToLoad);

