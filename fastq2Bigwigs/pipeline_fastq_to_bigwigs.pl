#!/usr/bin/env perl
# This script loads the treatment of data in the context of the lab core facility.
# It uses the input file files_.conf of format: FASTQFILESFOLDER;CHIPSEQ;GENOME;SPECIES;SINGLEEND;FRAGMENTSIZE;WIGFIXED;WIGVARIABLE;SPIKEIN 
#Nicolas Descostes august 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 12)
{
    die("Missing arguments for pipeline_fastq_to_bigwigs.pl: FASTQFILESFOLDER;CHIPSEQ;GENOME;SPECIES;SINGLEEND;FRAGMENTSIZE;WIGFIXED;WIGVARIABLE;SPIKEIN;MINCPU;MAXCPU\n\n");
}

my $fastqFilesFolder = $arguments_tab[0];
my $chipseq = $arguments_tab[1];
my $genome = $arguments_tab[2];
my $species = $arguments_tab[3];
my $singleEnd = $arguments_tab[4];
my $fragmentSize = $arguments_tab[5];
my $wigFixed = $arguments_tab[6];
my $wigVariable = $arguments_tab[7];
my $spikein = $arguments_tab[8];
my $analysis_name = $arguments_tab[9];
my $min_cpu_number = $arguments_tab[10];
my $max_cpu_number = $arguments_tab[11];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part1.R --fastqFilesFolder $fastqFilesFolder --chipseq $chipseq --genome $genome --species $species --singleEnd $singleEnd --fragmentSize $fragmentSize --wigFixed $wigFixed --wigVariable $wigVariable --spikein $spikein --analysisName $analysis_name --minCpuNumber $min_cpu_number --maxCpuNumber $max_cpu_number\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part1.R --fastqFilesFolder $fastqFilesFolder --chipseq $chipseq --genome $genome --species $species --singleEnd $singleEnd --fragmentSize $fragmentSize --wigFixed $wigFixed --wigVariable $wigVariable --spikein $spikein --analysisName $analysis_name --minCpuNumber $min_cpu_number --maxCpuNumber $max_cpu_number";
system($commandToLoad);




