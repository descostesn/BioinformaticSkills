#!/usr/bin/env perl
# This script performs gene set enrichment analysis from chip-seq peaks using the package Chip-enrich.
# It uses the input file files.conf of format: peaksfile;outputprefix;outputfolderbase;genomeversion;genesetsvector;locus;readlength;numpeakthresholdvec;fdr;mappability
#Nicolas Descostes dec 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 10)
{
	die("Missing arguments for generate_object.pl\n\n File should contain: peaksfile;outputprefix;outputfolderbase;genomeversion;genesetsvector;locus;readlength;numpeakthresholdvec;fdr;mappability\n\n");
}


my $peaksFile = $arguments_tab[0];
my $outputPrefix = $arguments_tab[1];
my $outputFolderBase = $arguments_tab[2];
my $genomeVersion = $arguments_tab[3];
my $genesSetsVector = $arguments_tab[4];
my $locus = $arguments_tab[5];
my $readLength = $arguments_tab[6];
my $numPeakThresholdVec = $arguments_tab[7];
my $FDR = $arguments_tab[8];
my $mappability = $arguments_tab[9];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/chipenrich_analysis.R --peaksFile $peaksFile --outputPrefix $outputPrefix --outputFolderBase $outputFolderBase --genomeVersion $genomeVersion --genesSetsVector $genesSetsVector --locus $locus --readLength $readLength --numPeakThresholdVec $numPeakThresholdVec --nbCpu $nbCPU --FDR $FDR --mappabilityOption $mappability\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/chipenrich_analysis.R --peaksFile $peaksFile --outputPrefix $outputPrefix --outputFolderBase $outputFolderBase --genomeVersion $genomeVersion --genesSetsVector $genesSetsVector --locus $locus --readLength $readLength --numPeakThresholdVec $numPeakThresholdVec --nbCpu $nbCPU --FDR $FDR --mappabilityOption $mappability";
system($commandToLoad);







