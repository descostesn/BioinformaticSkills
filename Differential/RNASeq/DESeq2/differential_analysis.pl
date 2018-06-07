#!/usr/bin/env perl
#This script performs a differential analysis
# It uses the input file_X.conf of format: BAMFILE;SINGLE;SAMPLEINFOFILE;NAMEOFREF;OUTPUTFOLDER;COMPARISONNAME;REFSEQANNO;EXPNAMES + SUROGATENB
#Nicolas Descostes February 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nb_cpu = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 9 && scalar(@arguments_tab) != 10)
{
	die("Missing arguments for differential_analysis.pl\n\n File should contain: BAMFILE;SINGLE;SAMPLEINFOFILE;NAMEOFREF;OUTPUTFOLDER;COMPARISONNAME;REFSEQANNO;EXPNAMES;OUTPUTFORMAT with or without SUROGATENB/performBatchCorrection (not both)\n\n");
}

my $bam_files_vec = $arguments_tab[0];
my $single = $arguments_tab[1];
my $samples_info_file = $arguments_tab[2];
my $name_of_reference = $arguments_tab[3];
my $output_folder = $arguments_tab[4];
my $comparison_name = $arguments_tab[5];
my $refseq_anno = $arguments_tab[6];
my $expnames_vec = $arguments_tab[7];
my $output_format = $arguments_tab[8];

print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 9)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format\n\n";
my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format";
system($commandToLoad);
	
}else{
	my $optional_argument = $arguments_tab[9];
	chomp($optional_argument);
	
	if ( $optional_argument =~ /^[0-9]+/ )
	{
		print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format --nbSurrogateVariable $optional_argument\n\n";
    	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format --nbSurrogateVariable $optional_argument";
		system($commandToLoad);
	}
	else
	{
		print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format --performBatchCorrection $optional_argument\n\n";
    	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/differential_analysis.R --bamFilesVec $bam_files_vec --nbCPU $nb_cpu --singleEnd $single --samplesInfoFile $samples_info_file --nameOfReference $name_of_reference --outputFolder $output_folder --comparisonName $comparison_name --refseqAnno $refseq_anno --expnamesVec $expnames_vec --outputFormat $output_format --performBatchCorrection $optional_argument";
		system($commandToLoad);
	}
}





