#!/usr/bin/env perl
#This script converts a wig file in variable steps to one in fixed steps
# It uses the input file files_wig_X.conf of format: WIGFOLDER;WIGFILE;NAMEWIG;BINNING;REFFILE
#Nicolas Descostes August 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
	die("Missing arguments for wigVariable2Fixed.pl\n\n File should contain: WIGFOLDER;WIGFILE;NAMEWIG;BINNING;REFFILE\n\n");
}

my $wigFolder = $arguments_tab[0];
my $wigFile = $arguments_tab[1];
my $nameWig = $arguments_tab[2];
my $binning = $arguments_tab[3];
my $refFile = $arguments_tab[4];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/wigVariable2Fixed.R --wigFolder $wigFolder --wigFile $wigFile --nameWig $nameWig --binning $binning --refFile $refFile\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/wigVariable2Fixed.R --wigFolder $wigFolder --wigFile $wigFile --nameWig $nameWig --binning $binning --refFile $refFile";
system($commandToLoad);




