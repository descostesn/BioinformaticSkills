#!/usr/bin/env perl
# This script load the pipeline by retrieving all parameters in a file files_pasha.conf
# It uses the input file files_pasha.conf of format: EXPNAME;FILEPATH;INPUTTYPE;ARTEFTH;BIN;REFFILE;PAIREDEND;MIDPT;SUBFOLD;REPORTFOLD;WIGFS;WIGVS;GFF;ELONGSIZE;NBCPU;KEEPTMP;ANNOGFF;MULTIFILE
#Nicolas Descostes april-june 2015

use strict;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);


if(scalar(@arguments_tab) != 18 && scalar(@arguments_tab) != 19)
{
	die("Missing arguments for pipeline_loader.pl\n\n\t File should contain: EXPNAME;FILEPATH;INPUTTYPE;ARTEFTH;BIN;REFFILE;PAIREDEND;MIDPT;SUBFOLD;REPORTFOLD;WIGFS;WIGVS;GFF;ELONGSIZE;NBCPU;KEEPTMP;ANNOGFF;MULTIFILE OR EXPNAME;FILEPATH;INPUTTYPE;ARTEFTH;BIN;REFFILE;PAIREDEND;MIDPT;SUBFOLD;REPORTFOLD;WIGFS;WIGVS;GFF;ELONGSIZE;NBCPU;KEEPTMP;INPUTCHRPREFIX;ANNOGFF;MULTIFILE");
}


print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 18)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/pipeline_loader.R --expName $arguments_tab[0] --inputFile $arguments_tab[1] --inputType $arguments_tab[2] --threshold $arguments_tab[3] --bins $arguments_tab[4] --annotationGenomeFiles $arguments_tab[5] --pairedEnds $arguments_tab[6] --midPoint $arguments_tab[7] --resultSubFolder $arguments_tab[8] --reportFilesSubFolder $arguments_tab[9] --WIGfs $arguments_tab[10] --WIGvs $arguments_tab[11] --GFF $arguments_tab[12] --elongationSize $arguments_tab[13] --ignoreChr \"random|hap\" --ignoreInsertsOver 500 --CPU $arguments_tab[14] --keepTemp $arguments_tab[15] --inputChrPrefix chr --inputChrSuffix '' --annotationFilesGFF $arguments_tab[16] --multireadSignalFile $arguments_tab[17]\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/pipeline_loader.R --expName $arguments_tab[0] --inputFile $arguments_tab[1] --inputType $arguments_tab[2] --threshold $arguments_tab[3] --bins $arguments_tab[4] --annotationGenomeFiles $arguments_tab[5] --pairedEnds $arguments_tab[6] --midPoint $arguments_tab[7] --resultSubFolder $arguments_tab[8] --reportFilesSubFolder $arguments_tab[9] --WIGfs $arguments_tab[10] --WIGvs $arguments_tab[11] --GFF $arguments_tab[12] --elongationSize $arguments_tab[13] --ignoreChr \"random|hap\" --ignoreInsertsOver 500 --CPU $arguments_tab[14] --keepTemp $arguments_tab[15] --inputChrPrefix chr --inputChrSuffix '' --annotationFilesGFF $arguments_tab[16] --multireadSignalFile $arguments_tab[17]";
	system($commandToLoad);
}
else
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/pipeline_loader.R --expName $arguments_tab[0] --inputFile $arguments_tab[1] --inputType $arguments_tab[2] --threshold $arguments_tab[3] --bins $arguments_tab[4] --annotationGenomeFiles $arguments_tab[5] --pairedEnds $arguments_tab[6] --midPoint $arguments_tab[7] --resultSubFolder $arguments_tab[8] --reportFilesSubFolder $arguments_tab[9] --WIGfs $arguments_tab[10] --WIGvs $arguments_tab[11] --GFF $arguments_tab[12] --elongationSize $arguments_tab[13] --ignoreChr \"random|hap\" --ignoreInsertsOver 500 --CPU $arguments_tab[14] --keepTemp $arguments_tab[15] --inputChrPrefix $arguments_tab[16] --inputChrSuffix '' --annotationFilesGFF $arguments_tab[17] --multireadSignalFile $arguments_tab[18]\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/pipeline_loader.R --expName $arguments_tab[0] --inputFile $arguments_tab[1] --inputType $arguments_tab[2] --threshold $arguments_tab[3] --bins $arguments_tab[4] --annotationGenomeFiles $arguments_tab[5] --pairedEnds $arguments_tab[6] --midPoint $arguments_tab[7] --resultSubFolder $arguments_tab[8] --reportFilesSubFolder $arguments_tab[9] --WIGfs $arguments_tab[10] --WIGvs $arguments_tab[11] --GFF $arguments_tab[12] --elongationSize $arguments_tab[13] --ignoreChr \"random|hap\" --ignoreInsertsOver 500 --CPU $arguments_tab[14] --keepTemp $arguments_tab[15] --inputChrPrefix $arguments_tab[16] --inputChrSuffix '' --annotationFilesGFF $arguments_tab[17] --multireadSignalFile $arguments_tab[18]";
	system($commandToLoad);
}

