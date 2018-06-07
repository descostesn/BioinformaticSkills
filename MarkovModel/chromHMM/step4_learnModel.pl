#!/usr/bin/env perl
#After having computed the binarized bed files, this script learn the HMM and output the different states.
# It uses the input file files_X.conf of format: OUTFILEID;GENOMEVERSION;INPUTDIR;OUTPUTDIR;NUMSTATES;SEGMENTATIONSIZE;CONVERGEDELTA;INFORMATIONSMOOTH;TYPEINITIALIZATION;MAXITERATIONS;STATEORDERING
#Nicolas Descostes DEC 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 11)
{
	die("Missing arguments for step4_learnModel.pl\n\n File should contain: OUTFILEID;GENOMEVERSION;INPUTDIR;OUTPUTDIR;NUMSTATES;SEGMENTATIONSIZE;CONVERGEDELTA;INFORMATIONSMOOTH;TYPEINITIALIZATION;MAXITERATIONS;STATEORDERING\n\n");
}


my $outfileID = $arguments_tab[0];
my $genomeVersion = $arguments_tab[1];
my $inputdir = $arguments_tab[2];
my $outputdir = $arguments_tab[3];
my $numstates = $arguments_tab[4];
my $segmentationSize = $arguments_tab[5];
my $convergedelta = $arguments_tab[6];
my $informationsmooth = $arguments_tab[7];
my $type_initialization = $arguments_tab[8];
my $maxiterations = $arguments_tab[9];
my $stateordering = $arguments_tab[10];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/step4_learnModel.R --outfileID $outfileID --genomeVersion $genomeVersion --nb_cpu $nbCPU --inputdir $inputdir --outputdir $outputdir --numstates $numstates --segmentationSize $segmentationSize --convergedelta $convergedelta --informationsmooth $informationsmooth --type_initialization $type_initialization --maxiterations $maxiterations --stateordering $stateordering\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/step4_learnModel.R --outfileID $outfileID --genomeVersion $genomeVersion --nb_cpu $nbCPU --inputdir $inputdir --outputdir $outputdir --numstates $numstates --segmentationSize $segmentationSize --convergedelta $convergedelta --informationsmooth $informationsmooth --type_initialization $type_initialization --maxiterations $maxiterations --stateordering $stateordering";
system($commandToLoad);




