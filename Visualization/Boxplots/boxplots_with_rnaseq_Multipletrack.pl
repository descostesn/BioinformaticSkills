#!/usr/bin/env perl
#This script generates boxplots for rna-seq experiments.
# It uses the input file filesX.conf of format: bamFile;expNames;stringIndexForGroups;outputFolder;refseqannofile;gffFile;singleEnd;filterGFFUniquely;indicateMeanWithSD;indicateMean;indicateMedian;indicateBoxplot;indicateJitter
#Nicolas Descostes nov 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 14)
{
	die("Missing arguments for boxplots_with_rnaseq.pl\n\n File should contain: bamFile;expNames;stringIndexForGroups;outputFolder;refseqannofile;gffFile;singleEnd;filterGFFUniquely;indicateMeanWithSD;indicateMean;indicateMedian;indicateBoxplot;indicateJitter;outputformat\n\n");
}

my $bam_file_vec = $arguments_tab[0];
my $expnames_vec = $arguments_tab[1];
my $string_index_for_groups = $arguments_tab[2];
my $output_folder = $arguments_tab[3];
my $refseq_anno_file = $arguments_tab[4];
my $gff_file = $arguments_tab[5];
my $single_end = $arguments_tab[6];
my $filter_gff_uniquely = $arguments_tab[7];
my $indicate_mean_with_sd = $arguments_tab[8];
my $indicate_mean = $arguments_tab[9];
my $indicate_median = $arguments_tab[10];
my $indicate_boxplot = $arguments_tab[11];
my $indicate_jitter = $arguments_tab[12];
my $output_format = $arguments_tab[13];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/boxplots_with_rnaseq_Multipletrack.R --bamFileVec $bam_file_vec --expnamesVec $expnames_vec --stringIndexForGroups $string_index_for_groups --outputFolder $output_folder --refseqAnnoFile $refseq_anno_file --gffFile $gff_file --singleEnd $single_end --filterGffUniquely $filter_gff_uniquely --indicateMeanWithSd $indicate_mean_with_sd --indicateMean $indicate_mean --indicateMedian $indicate_median --indicateBoxplot $indicate_boxplot --indicateJitter $indicate_jitter --nbCpu $nbCPU --outputFormat $output_format\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/boxplots_with_rnaseq_Multipletrack.R --bamFileVec $bam_file_vec --expnamesVec $expnames_vec --stringIndexForGroups $string_index_for_groups --outputFolder $output_folder --refseqAnnoFile $refseq_anno_file --gffFile $gff_file --singleEnd $single_end --filterGffUniquely $filter_gff_uniquely --indicateMeanWithSd $indicate_mean_with_sd --indicateMean $indicate_mean --indicateMedian $indicate_median --indicateBoxplot $indicate_boxplot --indicateJitter $indicate_jitter --nbCpu $nbCPU --outputFormat $output_format";
system($commandToLoad);








