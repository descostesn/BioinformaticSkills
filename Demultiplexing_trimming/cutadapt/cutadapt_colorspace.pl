#!/usr/bin/env perl
#This script removes adapters from color space sequences
# It uses the input file files_colorspace_X.conf of format: ADAPTERSEQ3PRIME;ADAPTERSEQ5PRIME;CSFASTAFILE;CSQUALFILE;OUTPUTFASTQ
#Nicolas Descostes July 2015

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
	die("Missing arguments for cutadapt_colorspace.pl\n\n File should contain: ADAPTERSEQ3PRIME;ADAPTERSEQ5PRIME;CSFASTAFILE;CSQUALFILE;OUTPUTFASTQ\n\n");
}


my $adapter_sequence3prime = $arguments_tab[0];
my $adapter_sequence5prime = $arguments_tab[1];
my $csfasta_file = $arguments_tab[2];
my $qual_file = $arguments_tab[3];
my $output_fastq = $arguments_tab[4];

		
print "This is job number $file_lineNumber \n";

print "cutadapt -c -a $adapter_sequence3prime -g $adapter_sequence5prime $csfasta_file $qual_file > $output_fastq\n";


my $commandToLoad = "cutadapt -c -a $adapter_sequence3prime -g $adapter_sequence5prime $csfasta_file $qual_file > $output_fastq";
system($commandToLoad);
















