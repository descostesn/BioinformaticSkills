#!/usr/bin/env perl
##########
# scripts that demultiplex and remove primers from sequencing files only if using Illumina primers.
# Nicolas Descostes, 01/27/2017
##########


use strict;
use warnings;

my $seq_file = $ARGV[0];

my $Illu1_primer = "CGATGT";
my $Illu2_primer = "TGACCA";
my $Illu3_primer = "ACAGTG";
my $Illu4_primer = "GCCAAT";
my $Illu5_primer = "CAGATC";
my $Illu6_primer = "CTTGTA";
my $Illu7_primer = "ATCACG";
my $Illu8_primer = "TTAGGC";
my $Illu9_primer = "ACTTGA";
my $Illu10_primer = "GATCAG";
my $Illu11_primer = "TAGCTT";
my $Illu12_primer = "GGCTAC";



if($seq_file !~ /fastq/)
{
    die("The input file should be in fastq format\n\n");
}

my @out_baseName_array = split(".fastq", $seq_file);
my $out_baseName = $out_baseName_array[0];

my $Illu1_file = join('', $out_baseName, "-Illu1.fastq");
my $Illu2_file = join('', $out_baseName, "-Illu2.fastq");
my $Illu3_file = join('', $out_baseName, "-Illu3.fastq");
my $Illu4_file = join('', $out_baseName, "-Illu4.fastq");
my $Illu5_file = join('', $out_baseName, "-Illu5.fastq");
my $Illu6_file = join('', $out_baseName, "-Illu6.fastq");
my $Illu7_file = join('', $out_baseName, "-Illu7.fastq");
my $Illu8_file = join('', $out_baseName, "-Illu8.fastq");
my $Illu9_file = join('', $out_baseName, "-Illu9.fastq");
my $Illu10_file = join('', $out_baseName, "-Illu10.fastq");
my $Illu11_file = join('', $out_baseName, "-Illu11.fastq");
my $Illu12_file = join('', $out_baseName, "-Illu12.fastq");
my $noprimer_file = join('', $out_baseName, "-noprimers.fastq");

open(FIC2,"> $Illu1_file") or die("\n\n Impossible to write the file $Illu1_file: $!\n\n");
open(FIC3,"> $Illu2_file") or die("\n\n Impossible to write the file $Illu2_file: $!\n\n");
open(FIC4,"> $Illu3_file") or die("\n\n Impossible to write the file $Illu3_file: $!\n\n");
open(FIC5,"> $Illu4_file") or die("\n\n Impossible to write the file $Illu4_file: $!\n\n");
open(FIC6,"> $Illu5_file") or die("\n\n Impossible to write the file $Illu5_file: $!\n\n");
open(FIC7,"> $Illu6_file") or die("\n\n Impossible to write the file $Illu6_file: $!\n\n");
open(FIC8,"> $Illu7_file") or die("\n\n Impossible to write the file $Illu7_file: $!\n\n");
open(FIC9,"> $Illu8_file") or die("\n\n Impossible to write the file $Illu8_file: $!\n\n");
open(FIC10,"> $Illu9_file") or die("\n\n Impossible to write the file $Illu9_file: $!\n\n");
open(FIC11,"> $Illu10_file") or die("\n\n Impossible to write the file $Illu10_file: $!\n\n");
open(FIC12,"> $Illu11_file") or die("\n\n Impossible to write the file $Illu11_file: $!\n\n");
open(FIC13,"> $Illu12_file") or die("\n\n Impossible to write the file $Illu12_file: $!\n\n");
open(FIC14,"> $noprimer_file") or die("\n\n Impossible to write the file $noprimer_file: $!\n\n");

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
        
        if($lines_block[1] =~ /^$Illu1_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu1_primer), $string_length-length($Illu1_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu1_primer), $string_length-length($Illu1_primer));         
            print(FIC2 join("\n", @lines_block));
            print(FIC2 "\n");   
        }
        elsif($lines_block[1] =~ /^$Illu2_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu2_primer), $string_length-length($Illu2_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu2_primer), $string_length-length($Illu2_primer));         
            print(FIC3 join("\n", @lines_block));
            print(FIC3 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu3_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu3_primer), $string_length-length($Illu3_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu3_primer), $string_length-length($Illu3_primer));         
            print(FIC4 join("\n", @lines_block));
            print(FIC4 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu4_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu4_primer), $string_length-length($Illu4_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu4_primer), $string_length-length($Illu4_primer));         
            print(FIC5 join("\n", @lines_block));
            print(FIC5 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu5_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu5_primer), $string_length-length($Illu5_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu5_primer), $string_length-length($Illu5_primer));         
            print(FIC6 join("\n", @lines_block));
            print(FIC6 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu6_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu6_primer), $string_length-length($Illu6_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu6_primer), $string_length-length($Illu6_primer));         
            print(FIC7 join("\n", @lines_block));
            print(FIC7 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu7_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu7_primer), $string_length-length($Illu7_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu7_primer), $string_length-length($Illu7_primer));         
            print(FIC8 join("\n", @lines_block));
            print(FIC8 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu8_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu8_primer), $string_length-length($Illu8_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu8_primer), $string_length-length($Illu8_primer));         
            print(FIC9 join("\n", @lines_block));
            print(FIC9 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu9_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu9_primer), $string_length-length($Illu9_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu9_primer), $string_length-length($Illu9_primer));         
            print(FIC10 join("\n", @lines_block));
            print(FIC10 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu10_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu10_primer), $string_length-length($Illu10_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu10_primer), $string_length-length($Illu10_primer));           
            print(FIC11 join("\n", @lines_block));
            print(FIC11 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu11_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu11_primer), $string_length-length($Illu11_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu11_primer), $string_length-length($Illu11_primer));           
            print(FIC12 join("\n", @lines_block));
            print(FIC12 "\n");
        }
        elsif($lines_block[1] =~ /^$Illu12_primer/)
        {
            my $string_length = length($lines_block[1]);
            $lines_block[1] = substr($lines_block[1], length($Illu12_primer), $string_length-length($Illu12_primer));
            $lines_block[3] = substr($lines_block[3], length($Illu12_primer), $string_length-length($Illu12_primer));           
            print(FIC13 join("\n", @lines_block));
            print(FIC13 "\n");
        }
        else
        {
            print(FIC14 join("\n", @lines_block));
            print(FIC14 "\n");
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
close(FIC9);
close(FIC10);
close(FIC11);
close(FIC12);
close(FIC13);
close(FIC14);
