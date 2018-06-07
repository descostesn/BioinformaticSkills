#!/usr/bin/env perl
#This script performs differential chip-seq peak analysis with MACS2 software
# It uses the input file files_X.conf of format: TREATMENTFILE;CONTROLFILE;EXPNAME;FORMAT;GENOMESIZE;TAGSIZE;QVALUE;MFOLD;ISBROAD
#Nicolas Descostes August 2016


######   !!!!!!!!!!!!!!!!! to continue: modify the peak detection to use bdg files !!!!!!!!!!!!!!!!!!!!!!!

use strict;
my $file_lineNumber = $ARGV[0];
my $input_conf_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_conf_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 9)
{
	die("Missing arguments for macs2diff.pl\n\n File should contain: TREATMENTFILE;CONTROLFILE;EXPNAME;FORMAT;GENOMESIZE;TAGSIZE;QVALUE;MFOLD;ISBROAD\n\n");
}

my $treatment_file = $arguments_tab[0]; # Bam file of the experiment
my $control_file = $arguments_tab[1]; # Bam file of the control input
my $experiment_name = $arguments_tab[2]; # Will be used for the output file name
my $format = $arguments_tab[3]; # "AUTO", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM" or "BOWTIE".
my $genome_size = $arguments_tab[4]; # Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
my $tag_size = $arguments_tab[5];    # Tag size. This will overide the auto detected tag size. DEFAULT: Not set
my $q_value = $arguments_tab[6];    # Minimum FDR (q-value) cutoff for peak detection. DEFAULT: 0.01
my $m_fold = $arguments_tab[7];    #Select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. The regions must be lower than upper limit, and higher than the lower limit. DEFAULT:10,30
my $is_broad = $arguments_tab[8];  # Should be set to 0 or 1


#Prediction of the extension size


if($is_broad)
{
	print "macs2 pedictd -t $treatment_file -c $control_file -n $experiment_name -f $format -g $genome_size -s $tag_size -q $q_value -m $m_fold --broad\n\n";
	my $commandToLoad = "macs2 pedictd -t $treatment_file -c $control_file -n $experiment_name -f $format -g $genome_size -s $tag_size -q $q_value -m $m_fold --broad";
	system($commandToLoad);
	
}else{
	
	print "macs2 pedictd -t $treatment_file -c $control_file -n $experiment_name -f $format -g $genome_size -s $tag_size -q $q_value -m $m_fold\n\n";
	my $commandToLoad = "macs2 pedictd -t $treatment_file -c $control_file -n $experiment_name -f $format -g $genome_size -s $tag_size -q $q_value -m $m_fold";
	system($commandToLoad);
}



