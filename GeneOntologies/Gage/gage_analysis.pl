#!/usr/bin/env perl
#This script performs GSEA analysis with GAGE
# It uses the input file_X.conf of format: SPECIESNAME;DOWNREGULATEDFILE;UPREGULATEDFILE;OUTPUTFILENAME;INPUTTYPE;NBCol with options MSigDBFile;MSigDBNames
#Nicolas Descostes nov 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 6 && scalar(@arguments_tab) != 8)
{
	die("Missing arguments for gage_analysis.pl\n\n File should contain: SPECIESNAME;DOWNREGULATEDFILE;UPREGULATEDFILE;OUTPUTFILENAME;INPUTTYPE;NBCol with options MSigDBFile;MSigDBNames\n\n");
}


my $species_name = $arguments_tab[0];
my $downregulated_file = $arguments_tab[1];
my $upregulated_file = $arguments_tab[2];
my $output_file_name = $arguments_tab[3];
my $input_type = $arguments_tab[4];
my $nb_cols = $arguments_tab[5];

print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 6)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gage_analysis.R --speciesName $species_name --downregulatedFile $downregulated_file --upregulatedFile $upregulated_file --outputFileName $output_file_name --inputType $input_type --nbCols $nb_cols\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gage_analysis.R --speciesName $species_name --downregulatedFile $downregulated_file --upregulatedFile $upregulated_file --outputFileName $output_file_name --inputType $input_type --nbCols $nb_cols";
	system($commandToLoad);
	
}else{

	my $misigDBFileVec  = $arguments_tab[6];
	my $msigDBNameVec  = $arguments_tab[7];
	
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gage_analysis.R --speciesName $species_name --downregulatedFile $downregulated_file --upregulatedFile $upregulated_file --outputFileName $output_file_name --inputType $input_type --nbCols $nb_cols --msigDBFileVec $misigDBFileVec --msigDBNameVec $msigDBNameVec\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/gage_analysis.R --speciesName $species_name --downregulatedFile $downregulated_file --upregulatedFile $upregulated_file --outputFileName $output_file_name --inputType $input_type --nbCols $nb_cols --msigDBFileVec $misigDBFileVec --msigDBNameVec $msigDBNameVec";
	system($commandToLoad);
}
