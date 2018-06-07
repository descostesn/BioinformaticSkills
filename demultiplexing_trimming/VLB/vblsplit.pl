#!/usr/bin/env perl
##########
# scripts that demultiplex and remove primers from sequencing files only if using VBL primers.
# Nicolas Descostes, 06/09/2015
##########


use strict;
use warnings;

my $seq_file = $ARGV[0];
my $VLB1_primer = "AT";
my $VLB2_primer = "CAT";
my $VLB3_primer = "GCAT";
my $VLB4_primer = "TGCAT";
my $VLB5_primer = "CTGCAT";
my $VLB6_primer = "GCTGCAT";

if($seq_file !~ /fastq/)
{
	die("The input file should be in fastq format\n\n");
}

my @out_baseName_array = split(".fastq", $seq_file);
my $out_baseName = $out_baseName_array[0];

my $vlb1_file = join('', $out_baseName, "-VLB1.fastq");
my $vlb2_file = join('', $out_baseName, "-VLB2.fastq");
my $vlb3_file = join('', $out_baseName, "-VLB3.fastq");
my $vlb4_file = join('', $out_baseName, "-VLB4.fastq");
my $vlb5_file = join('', $out_baseName, "-VLB5.fastq");
my $vlb6_file = join('', $out_baseName, "-VLB6.fastq");
my $noprimer_file = join('', $out_baseName, "-noprimers.fastq");

open(FIC2,"> $vlb1_file") or die("\n\n Impossible to write the file $vlb1_file: $!\n\n");
open(FIC3,"> $vlb2_file") or die("\n\n Impossible to write the file $vlb2_file: $!\n\n");
open(FIC4,"> $vlb3_file") or die("\n\n Impossible to write the file $vlb3_file: $!\n\n");
open(FIC5,"> $vlb4_file") or die("\n\n Impossible to write the file $vlb4_file: $!\n\n");
open(FIC6,"> $vlb5_file") or die("\n\n Impossible to write the file $vlb5_file: $!\n\n");
open(FIC7,"> $vlb6_file") or die("\n\n Impossible to write the file $vlb6_file: $!\n\n");
open(FIC8,"> $noprimer_file") or die("\n\n Impossible to write the file $noprimer_file: $!\n\n");

open(FIC1,"< $seq_file") or die("\n\n Impossible to open the file $seq_file: $!\n\n");

my $countline = 0;
my @lines_block = ();
my $line_to_write = "";

while(defined(my $line =<FIC1>))
{
	$countline++;
	chomp($line);
	push @lines_block, $line; 
	
	if($countline == 4)
	{
		
		if($lines_block[1] =~ /^$VLB1_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB1_primer), $string_length-length($VLB1_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB1_primer), $string_length-length($VLB1_primer));			
			print(FIC2 join("\n", @lines_block));
			print(FIC2 "\n");	
		}
		elsif($lines_block[1] =~ /^$VLB2_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB2_primer), $string_length-length($VLB2_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB2_primer), $string_length-length($VLB2_primer));			
			print(FIC3 join("\n", @lines_block));
			print(FIC3 "\n");
		}
		elsif($lines_block[1] =~ /^$VLB3_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB3_primer), $string_length-length($VLB3_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB3_primer), $string_length-length($VLB3_primer));			
			print(FIC4 join("\n", @lines_block));
			print(FIC4 "\n");
		}
		elsif($lines_block[1] =~ /^$VLB4_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB4_primer), $string_length-length($VLB4_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB4_primer), $string_length-length($VLB4_primer));			
			print(FIC5 join("\n", @lines_block));
			print(FIC5 "\n");
		}
		elsif($lines_block[1] =~ /^$VLB5_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB5_primer), $string_length-length($VLB5_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB5_primer), $string_length-length($VLB5_primer));			
			print(FIC6 join("\n", @lines_block));
			print(FIC6 "\n");
		}
		elsif($lines_block[1] =~ /^$VLB6_primer/)
		{
			my $string_length = length($lines_block[1]);
			$lines_block[1] = substr($lines_block[1], length($VLB6_primer), $string_length-length($VLB6_primer));
			$lines_block[3] = substr($lines_block[3], length($VLB6_primer), $string_length-length($VLB6_primer));			
			print(FIC7 join("\n", @lines_block));
			print(FIC7 "\n");
		}
		else
		{
			print(FIC8 join("\n", @lines_block));
			print(FIC8 "\n");
		}
		
		$countline = 0;
		@lines_block = ();
	}
}


close(FIC1);
close(FIC2);
close(FIC3);
close(FIC4);
close(FIC5);
close(FIC6);
close(FIC7);
close(FIC8);