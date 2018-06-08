#!/usr/bin/env perl
# This script extract coordinates for each pattern
# It uses the input file files_X.conf of format: FOLDERVEC
#Nicolas Descostes December 2015

use strict;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);


if(scalar(@arguments_tab) != 1)
{
	die("Missing arguments for extract_from_allpatterns.pl\n\n\t File should contain: FOLDERVEC");
}

my $working_dir = $arguments_tab[0];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/extract_from_allpatterns.R --workingDirVec $working_dir\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/extract_from_allpatterns.R --workingDirVec $working_dir";
system($commandToLoad); 


	








 

