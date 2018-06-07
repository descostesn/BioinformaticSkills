################
#This script enables to convert a bed file to a gff file (for the pipeline)
# command: perl bedTogff.pl path_to_input/ inputFile exp_Name
#exp_Name: the type of data (3rd column of the gff file)
#Descostes 15/09/2011
################

#!/usr/bin/perl -w
use strict;
use Benchmark;


my $pathToInput = $ARGV[0];
my $inputFile = $ARGV[1];
my $parameter_anno = $ARGV[2];

my @outputFileTab = split(/\./,$inputFile);
my $outputFile = $outputFileTab[0].".gff";

print "$outputFile\n";

if(!defined($inputFile) || !defined($pathToInput) || !defined($parameter_anno))
{
	die("\n\n invalid command: perl bedTogff.pl path_to_input/ inputFile exp_Name\n\n");
}



open(FIC1,"<$pathToInput$inputFile") or die("\n\n Impossible to open the file $inputFile: $!\n\n");
open(FIC2,">$pathToInput$outputFile") or die("\n\n Impossible to open the file $outputFile: $!\n\n");

my $count = 1;

while(defined(my $line =<FIC1>))
{
	chomp($line);
	my @informationBed = split(/\t/,$line);  
	my $chr = $informationBed[0];
	my $nameAnno = "";
	if (length($informationBed[3]) > 0){
		$nameAnno = $informationBed[3];
	}else{
		$nameAnno = "anno$count";
	    $count++;
	}
	my $typeAnno = $parameter_anno;
	my $start = $informationBed[1];
	my $end = $informationBed[2];
	my $score = 0;
	my $strand = "";
	if (length($informationBed[5]) > 0){
		$strand = $informationBed[5];
	}else{
		$strand = '+';
	}
	my $frame = ".";
	my $group = ".";
	
	print(FIC2 "$chr\t$typeAnno\t$nameAnno\t$start\t$end\t$score\t$strand\t$frame\t$group\n");

}


close(FIC1);
close(FIC2);