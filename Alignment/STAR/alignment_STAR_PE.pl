#!/usr/bin/env perl
# This script aligns fastq files with STAR aligner
# It uses the input file of format: GENOMEINDEXPATH;FASTQFILE1;FASTQFILE2;LIMITNBMULTIREADS;OVERHANGUNANNOTATEDJUNCTION;OVERHANGANNOTATEDJUNCTION;MININTRONLENGTH;MAXINTRONLENGTH;OUTPUTFILENAMEPREFIX;QUALITYTHRESHOLD;NORMWIGGLE;READLENGTH;STRANDWIG
#Nicolas Descostes March 2018

use strict;
use File::Basename;

my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];
my $nb_cpu = $ARGV[2];
my $genome_index = $ARGV[3];

#Retrieve the parameters
my $file_param = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_param;
my @param_tab = split(';', $file_param);

if(scalar(@param_tab) != 12)
{
    die("Missing arguments for alignment_STAR_PE.pl\n\n\t input file conf should have FASTQFILE1;FASTQFILE2;LIMITNBMULTIREADS;OVERHANGUNANNOTATEDJUNCTION;OVERHANGANNOTATEDJUNCTION;MININTRONLENGTH;MAXINTRONLENGTH;OUTPUTFILENAMEPREFIX;QUALITYTHRESHOLD;NORMWIGGLE;READLENGTH;STRANDWIG\n\n");    
}



my $fastq_file1 = $param_tab[0];
my $fastq_file2 = $param_tab[1];
my $limit_nb_multireads = $param_tab[2];
my $overhang_unannotated_junctions = $param_tab[3];     # Default value: 8 
my $overhang_annotated_junctions = $param_tab[4];       # Default value: 1
my $min_intron_length = $param_tab[5];                  # Default value: 20
my $max_intron_length = $param_tab[6];                  # Default value: 1000000
my $quality_threshold = $param_tab[8];                  # Default value: 255
my $norm_wiggle = $param_tab[9];                        # RPM or None
my $read_length = $param_tab[10];
my $strand_wig = $param_tab[11];

my $output_folder = join("", dirname($fastq_file1), "/");
my $output_file_namePrefix = join("", $output_folder, "STAR_", $param_tab[7], "/");
mkdir $output_file_namePrefix or die "Error creating directory: $output_file_namePrefix";
$output_file_namePrefix = join("", $output_file_namePrefix, $param_tab[7]);

print "This is job number $file_lineNumber \n";

my $outFilter_mismatchNMax = (0.06 * $read_length);

print "STAR --runThreadN $nb_cpu --genomeDir $genome_index --readFilesIn $fastq_file1 $fastq_file2 --outFilterMultimapNmax $limit_nb_multireads --alignSJoverhangMin $overhang_unannotated_junctions --alignSJDBoverhangMin $overhang_annotated_junctions --outFilterMismatchNmax $outFilter_mismatchNMax --alignIntronMin $min_intron_length --alignIntronMax $max_intron_length --outFileNamePrefix $output_file_namePrefix --outSAMmapqUnique $quality_threshold --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigNorm $norm_wiggle --outWigStrand $strand_wig\n";

my $commandToLoad = "STAR --runThreadN $nb_cpu --genomeDir $genome_index --readFilesIn $fastq_file1 $fastq_file2 --outFilterMultimapNmax $limit_nb_multireads --alignSJoverhangMin $overhang_unannotated_junctions --alignSJDBoverhangMin $overhang_annotated_junctions --outFilterMismatchNmax $outFilter_mismatchNMax --alignIntronMin $min_intron_length --alignIntronMax $max_intron_length --outFileNamePrefix $output_file_namePrefix --outSAMmapqUnique $quality_threshold --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigNorm $norm_wiggle --outWigStrand $strand_wig";
system($commandToLoad);
 




 



