#!/usr/bin/env perl
#This script detects q-quadruplexes in given sequences
# It uses the input file files_X.conf of format: FASTAFILE;SEQUENCELENGTH;OUTPUTFOLDER;EXPNAME;OBJECTVEC;SORTGFF;SPACESIZE;INTERVUP;INTERVDOWN;BINSIZE;PROLENGTH
#Nicolas Descostes May 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];
my $nbcpu = $ARGV[2];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 11)
{
	die("Missing arguments for q-quadruplex.pl\n\n File should contain: FASTAFILE;SEQUENCELENGTH;OUTPUTFOLDER;EXPNAME;OBJECTVEC;SORTGFF;SPACESIZE;INTERVUP;INTERVDOWN;BINSIZE;PROLENGTH\n\n");
}

my $fasta_file = $arguments_tab[0];
my $sequence_length = $arguments_tab[1];
my $outputFolder = $arguments_tab[2];
my $expname = $arguments_tab[3];
my $object_vec = $arguments_tab[4];
my $sort_gff = $arguments_tab[5];
my $spacesize = $arguments_tab[6];
my $intervalupstream = $arguments_tab[7];
my $intervaldownstream = $arguments_tab[8];
my $binsize = $arguments_tab[9];
my $profilelength = $arguments_tab[10];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/g-quadruplex.R --fastaFileVec $fasta_file --seqLength $sequence_length --outputFolder $outputFolder --expname $expname --nbCpu $nbcpu --objectVec $object_vec --sortedCoordGffVec $sort_gff --spaceSizeObject $spacesize --intervalUpstream $intervalupstream --intervalDownstream $intervaldownstream --binsize $binsize --profileLength $profilelength\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/g-quadruplex.R --fastaFileVec $fasta_file --seqLength $sequence_length --outputFolder $outputFolder --expname $expname --nbCpu $nbcpu --objectVec $object_vec --sortedCoordGffVec $sort_gff --spaceSizeObject $spacesize --intervalUpstream $intervalupstream --intervalDownstream $intervaldownstream --binsize $binsize --profileLength $profilelength";
system($commandToLoad);


