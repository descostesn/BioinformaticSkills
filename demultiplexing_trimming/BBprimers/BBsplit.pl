#!/usr/bin/env perl
##########
# scripts that demultiplex and remove primers from sequencing files only if using BB primers.
# Nicolas Descostes, october 2015
##########


use strict;
use warnings;

my $seq_file = $ARGV[0];
my $BB1_primer = "GA";
my $BB2_primer = "CCA";
my $BB3_primer = "AGCA";
my $BB4_primer = "TTGCA";
my $BB5_primer = "ACTGCA";
my $BB6_primer = "TGCTGCA";

if($seq_file !~ /fastq/)
{
    die("The input file should be in fastq format\n\n");
}

my @out_baseName_array = split(".fastq", $seq_file);
my $out_baseName = $out_baseName_array[0];

my $bb1_file = join('', $out_baseName, "-BB1.fastq");
my $bb2_file = join('', $out_baseName, "-BB2.fastq");
my $bb3_file = join('', $out_baseName, "-BB3.fastq");
my $bb4_file = join('', $out_baseName, "-BB4.fastq");
my $bb5_file = join('', $out_baseName, "-BB5.fastq");
my $bb6_file = join('', $out_baseName, "-BB6.fastq");
my $noprimer_file = join('', $out_baseName, "-noprimers.fastq");

open(FIC2,"> $bb1_file") or die("\n\n Impossible to write the file $bb1_file: $!\n\n");
open(FIC3,"> $bb2_file") or die("\n\n Impossible to write the file $bb2_file: $!\n\n");
open(FIC4,"> $bb3_file") or die("\n\n Impossible to write the file $bb3_file: $!\n\n");
open(FIC5,"> $bb4_file") or die("\n\n Impossible to write the file $bb4_file: $!\n\n");
open(FIC6,"> $bb5_file") or die("\n\n Impossible to write the file $bb5_file: $!\n\n");
open(FIC7,"> $bb6_file") or die("\n\n Impossible to write the file $bb6_file: $!\n\n");
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
        
        if($lines_block[1] =~ /^$BB1_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB1_primer), $string_length-length($BB1_primer));
            $lines_block[3] = substr($lines_block[3], length($BB1_primer), $string_length-length($BB1_primer));         
            print(FIC2 join("\n", @lines_block));
            print(FIC2 "\n");   
        }
        elsif($lines_block[1] =~ /^$BB2_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB2_primer), $string_length-length($BB2_primer));
            $lines_block[3] = substr($lines_block[3], length($BB2_primer), $string_length-length($BB2_primer));         
            print(FIC3 join("\n", @lines_block));
            print(FIC3 "\n");
        }
        elsif($lines_block[1] =~ /^$BB3_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB3_primer), $string_length-length($BB3_primer));
            $lines_block[3] = substr($lines_block[3], length($BB3_primer), $string_length-length($BB3_primer));         
            print(FIC4 join("\n", @lines_block));
            print(FIC4 "\n");
        }
        elsif($lines_block[1] =~ /^$BB4_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB4_primer), $string_length-length($BB4_primer));
            $lines_block[3] = substr($lines_block[3], length($BB4_primer), $string_length-length($BB4_primer));         
            print(FIC5 join("\n", @lines_block));
            print(FIC5 "\n");
        }
        elsif($lines_block[1] =~ /^$BB5_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB5_primer), $string_length-length($BB5_primer));
            $lines_block[3] = substr($lines_block[3], length($BB5_primer), $string_length-length($BB5_primer));         
            print(FIC6 join("\n", @lines_block));
            print(FIC6 "\n");
        }
        elsif($lines_block[1] =~ /^$BB6_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($BB6_primer), $string_length-length($BB6_primer));
            $lines_block[3] = substr($lines_block[3], length($BB6_primer), $string_length-length($BB6_primer));         
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