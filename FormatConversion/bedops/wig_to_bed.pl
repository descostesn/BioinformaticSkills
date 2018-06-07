#!/usr/bin/env perl
#This script converts the wig file to bed file
# It uses the input file X.conf of format: INPUTWIG;OUTPUTBED
#Nicolas Descostes march 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_param = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_param;
my @param_tab = split(';', $file_param);

if(scalar(@param_tab) != 2)
{
	die("Missing arguments for wig_to_bed.pl\n\n\t input file conf should have INPUTWIG;OUTPUTBED\n\n");	
}

my $input_wig = $param_tab[0];
my $output_bed = $param_tab[1];
						
print "This is job number $file_lineNumber \n";

print "wig2bed < $input_wig > $output_bed\n";

my $commandToLoad = "wig2bed < $input_wig > $output_bed";
system($commandToLoad);

 