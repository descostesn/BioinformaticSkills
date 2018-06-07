#!/usr/bin/env perl
#  Preprocess bam files into binary files used by PeakSeq by pipping with samtools.
# It uses the input file files_wig_X.conf of format: BAMFILE;EXPNAME;OUTPUTFOLDER 
#Nicolas Descostes march 2017

use strict;
use File::Basename qw( fileparse );
use File::Path qw( make_path );
use File::Spec;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 3)
{
	die("Missing arguments for preprocess_bam_files_PeakSeq.pl\n\n File should contain: BAMFILE;EXPNAME;OUTPUTFOLDER\n\n");
}

my $file_bam = $arguments_tab[0]; 
my $expname = $arguments_tab[1];
my $output_folder = $arguments_tab[2];

my $filenamePath = $output_folder . $expname;

if ( !-d $filenamePath ) {
    make_path $filenamePath or die "Failed to create path: $filenamePath";
}


print "This is job number $file_lineNumber \n";

print "samtools view $file_bam | /ifs/home/descon01/programs/PeakSeq/bin/PeakSeq -preprocess SAM stdin $filenamePath\n\n";

my $commandToLoad = "samtools view $file_bam | /ifs/home/descon01/programs/PeakSeq/bin/PeakSeq -preprocess SAM stdin $filenamePath";
system($commandToLoad);
