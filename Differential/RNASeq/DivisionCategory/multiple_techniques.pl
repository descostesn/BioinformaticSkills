#!/usr/bin/env perl
# This script plots the distribution of gene expression from a bam file. define different classes with different classification techniques, and output bed files for each category.
# It uses the input file files_X.conf of format: BAMFILE;GTFFILE;SINGLEEND;OUTPUTFOLDER;BEDANNO;NBGROUPS
#Nicolas Descostes feb 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 6)
{
	die("Missing arguments for multiple_techniques.pl: BAMFILE;GTFFILE;SINGLEEND;OUTPUTFOLDER;BEDANNO;NBGROUPS\n\n");
}

my $bam_file= $arguments_tab[0];
my $gtf_file= $arguments_tab[1];
my $single_end= $arguments_tab[2];
my $output_folder= $arguments_tab[3];
my $bed_anno= $arguments_tab[4];
my $nb_of_groups= $arguments_tab[5];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/multiple_techniques.R --bamFile $bam_file --gtfFile $gtf_file --singleEnd $single_end --outputFolder $output_folder --bedAnno $bed_anno --nbOfGroups $nb_of_groups\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/multiple_techniques.R --bamFile $bam_file --gtfFile $gtf_file --singleEnd $single_end --outputFolder $output_folder --bedAnno $bed_anno --nbOfGroups $nb_of_groups";
system($commandToLoad);










