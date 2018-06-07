###############
# This script binarizes bed files that were obtained from bam files in step1.
# conf file should contain: CONTROLDIRECTORY;GENOMEVERSION;INPUTBEDDIR;CELLMARKFILETABLE;OUTPUTFOLDER;SEGMENTATIONSIZE;FOLDTHRES;SIGNALTHRES;SHIFT;POISSONTHRES;PSEUDOCOUNTCTRL;FLANKWIDTHCONTROL
# Descostes DEC 2016
###############

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 12)
{
	die("Missing arguments for step3_binarizeBedFiles.pl\n\n File should contain: CONTROLDIRECTORY;GENOMEVERSION;INPUTBEDDIR;CELLMARKFILETABLE;OUTPUTFOLDER;SEGMENTATIONSIZE;FOLDTHRES;SIGNALTHRES;SHIFT;POISSONTHRES;PSEUDOCOUNTCTRL;FLANKWIDTHCONTROL\n\n");
}

my $controlDirectory = $arguments_tab[0];
my $genomeVersion = $arguments_tab[1];
my $inputbeddir = $arguments_tab[2];
my $cellmarkfiletable = $arguments_tab[3];
my $outputFolder = $arguments_tab[4];
my $segmentationSize = $arguments_tab[5];
my $foldthresh = $arguments_tab[6];
my $signalthresh = $arguments_tab[7];
my $shift = $arguments_tab[8];
my $poissonthresh = $arguments_tab[9];
my $pseudocountcontrol = $arguments_tab[10];
my $flankwidthcontrol = $arguments_tab[11];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/step3_binarizeBedFiles.R --controlDirectory $controlDirectory --genomeVersion $genomeVersion --inputbeddir $inputbeddir --cellmarkfiletable $cellmarkfiletable --outputFolder $outputFolder --segmentationSize $segmentationSize --foldthresh $foldthresh --signalthresh $signalthresh --shift $shift --poissonthresh $poissonthresh --pseudocountcontrol $pseudocountcontrol --flankwidthcontrol $flankwidthcontrol\n\n";
my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/step3_binarizeBedFiles.R --controlDirectory $controlDirectory --genomeVersion $genomeVersion --inputbeddir $inputbeddir --cellmarkfiletable $cellmarkfiletable --outputFolder $outputFolder --segmentationSize $segmentationSize --foldthresh $foldthresh --signalthresh $signalthresh --shift $shift --poissonthresh $poissonthresh --pseudocountcontrol $pseudocountcontrol --flankwidthcontrol $flankwidthcontrol";
system($commandToLoad);





