#!/usr/bin/env perl
# This script performs motif discovery with meme using a background
# It uses the input file files.conf of format: see script]
#Nicolas Descostes March 2018

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 8)
{
    die("Missing arguments for meme.pl\n\n File should contain: see script\n\n");
}


my $fasta_file  = $arguments_tab[0];              #file containing sequences in FASTA format
my $output_dir = $arguments_tab[1];               # name of directory for output files will replace existing directory
my $nmotifs = $arguments_tab[2];                  # maximum number of motifs to find
my $minw = $arguments_tab[3];                     # minumum motif width. default is 8
my $maxw = $arguments_tab[4];                     # maximum motif width default is 50
my $mod = $arguments_tab[5];
my $maxSize = $arguments_tab[6];
my $backgroundModel = $arguments_tab[7];


print "This is job number $file_lineNumber \n";

print "meme $fasta_file -oc $output_dir -dna -mod $mod -minw $minw -maxw $maxw -nmotifs $nmotifs -maxsize $maxSize -bfile $backgroundModel -p $nbCPU -V\n\n";

my $commandToLoad = "meme $fasta_file -oc $output_dir -dna -mod $mod -minw $minw -maxw $maxw -nmotifs $nmotifs -maxsize $maxSize -bfile $backgroundModel -p $nbCPU -V";
system($commandToLoad);

