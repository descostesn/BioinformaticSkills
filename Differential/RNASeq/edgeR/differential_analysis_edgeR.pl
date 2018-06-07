#!/usr/bin/env perl
# This script performs differential analysis using edgeR package.
# It uses the input file files_correl_X.conf of format: BAMFILEVEC;REFSEQANNO;SINGLEEND;EXPNAMESVEC;CONDITIONVEC;OUTPUTFOLDER;SPECIES;GFFFILE;GROUPNAMES;CELL 
#Nicolas Descostes october 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 10)
{
	die("Missing arguments for differential_analysis_edgeR.pl\n\n File should contain: BAMFILEVEC;REFSEQANNO;SINGLEEND;EXPNAMESVEC;CONDITIONVEC;OUTPUTFOLDER;SPECIES;GFFFILE;GROUPNAMES;CELL\n\n");
}

my $bam_files_vec  = $arguments_tab[0];
my $refseq_anno  = $arguments_tab[1];
my $single_end  = $arguments_tab[2];
my $expnames_vec  = $arguments_tab[3];
my $conditions_vec  = $arguments_tab[4];
my $output_folder  = $arguments_tab[5];
my $species  = $arguments_tab[6];
my $gffFile = $arguments_tab[7];
my $groupnames = $arguments_tab[8];
my $cell = $arguments_tab[9];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis_edgeR.R --bamFilesVec $bam_files_vec --refseqAnno $refseq_anno --singleEnd $single_end --expnamesVec $expnames_vec --conditionsVec $conditions_vec --outputFolder $output_folder --species $species --gffFile $gffFile --groupnames $groupnames --cell $cell\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis_edgeR.R --bamFilesVec $bam_files_vec --refseqAnno $refseq_anno --singleEnd $single_end --expnamesVec $expnames_vec --conditionsVec $conditions_vec --outputFolder $output_folder --species $species --gffFile $gffFile --groupnames $groupnames --cell $cell";
system($commandToLoad);












