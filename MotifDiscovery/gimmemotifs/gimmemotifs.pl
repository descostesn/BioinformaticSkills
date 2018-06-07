#!/usr/bin/env perl
# This script performs motif discovery with gimmemotifs
# It uses the input file files.conf of format: nameexp;sizeanalysis;genomeversion;fractionused;widthmotif;enrichmentcutoff;pvalue;toolsused;markovorder;fastafile
#Nicolas Descostes dec 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 10)
{
	die("Missing arguments for gimmemotifs.pl\n\n File should contain: nameexp;sizeanalysis;genomeversion;fractionused;widthmotif;enrichmentcutoff;pvalue;toolsused;markovorder;fastafile\n\n");
}


my $name_experiment = $arguments_tab[0];
my $size_analysis = $arguments_tab[1];
my $genome_version = $arguments_tab[2];
my $fraction_used = $arguments_tab[3];
my $width_motif = $arguments_tab[4];
my $enrichment_cutoff = $arguments_tab[5];
my $pvalue = $arguments_tab[6];
my $tools_used = $arguments_tab[7];
my $markov_order = $arguments_tab[8];
my $fasta_file = $arguments_tab[9];

print "This is job number $file_lineNumber \n";

print "gimme motifs --keepintermediate -n $name_experiment -a $size_analysis -g $genome_version -f $fraction_used -w $width_motif -e $enrichment_cutoff -p $pvalue -t $tools_used -m $markov_order $fasta_file\n\n";

my $commandToLoad = "gimme motifs --keepintermediate -n $name_experiment -a $size_analysis -g $genome_version -f $fraction_used -w $width_motif -e $enrichment_cutoff -p $pvalue -t $tools_used -m $markov_order $fasta_file";
system($commandToLoad);












