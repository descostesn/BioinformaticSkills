#!/usr/bin/env perl
#This script performs peak detection with macs software
# It uses the input file files_MACS_X.conf of format: BAMFILEVEC;INPUTFILEVEC;EXPNAME;OUTPUTFOLDER;FORMAT;GENOMESIZE;TAGSIZE;QVALUE;ELONGATIONSIZE;ARTEFACTTHRESHOLD
# default p-value is 1e-5
#Nicolas Descostes february 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $input_conf_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_conf_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 10)
{
	die("Missing arguments for macs_peakcalling.pl\n\n File should contain: BAMFILEVEC;INPUTFILEVEC;EXPNAME;OUTPUTFOLDER;FORMAT;GENOMESIZE;TAGSIZE;QVALUE;ELONGATIONSIZE;ARTEFACTTHRESHOLD\n\n");
}


#$bam_file_vector: can be one or more bam files
#$input_file_vector: corresponding input files to the bam files
#$experiment_name: name that will be used as prefix for the output files
#$output_folder: folder to which all output files will be written
#$format: format of the files passed to -t. If several format, used 'AUTO'. Accepted format are "ELAND", "BED", "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET" (for pair-end tags), "SAM", "BAM", "BOWTIE" or "BAMPE".
#$genome_size: Size of the used genome.
# $tag_size: Size of the reads used, can be found in the fastq files. 
# $sonication_size: the band width that is used to scan the genome, it corresponds to the sonication fragment size expected from wet experiment (check the estimated size from the pipeline).
# $qvalue: Qvalue for peak detection.
# $elongation_size: Elongation size by which tags should be extended. Determine by Pasha, see log reports.
# $artefact_threshold: Threshold to remove duplicated tags as it is used in Pasha (not number of tags but actual threshold). See log reports of pasha to determine this parameter. 
# !!! option -m not used: This parameter is used to select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. The regions must be lower than upper limit, and higher than the lower limit of fold enrichment. DEFAULT:5,50 means using all regions not too low (>5) and not too high (<50) to build paired-peaks model. If MACS can not find more than 100 regions to build model, it will use the --extsize parameter to continue the peak detection ONLY if --fix-bimodal is set.
# !!! option --nolambda not used: With this flag on, MACS will use the background lambda as local lambda. This means MACS will not consider the local bias at peak candidate regions.
# !!! option --fix-bimodal not used: Whether turn on the auto paired-peak model process. If it's set, when MACS failed to build paired model, it will use the nomodel settings, the '--extsize' parameter to extend each tags. If set, MACS will be terminated if paried-peak model is failed.
# !!! option --bw $sonication_size: not used if --nomodel option is declared. If the model of macs wants to be used, check the documentation at https://github.com/taoliu/MACS.
# !!! option --to-large not used: When set, linearly scale the smaller dataset to the same depth as larger dataset, by default, the larger dataset will be scaled towards the smaller dataset. Beware, to scale up small data would cause more false positives.
# !!! option --call-summits used: to find subpeaks, see documentation
#!!! option --slocal, --llocal not used: These two parameters control which two levels of regions will be checked around the peak regions to calculate the maximum lambda as local lambda. By default, MACS considers 1000bp for small local region(--slocal), and 10000bps for large local region(--llocal) which captures the bias from a long range effect like an open chromatin domain. You can tweak these according to your project. Remember that if the region is set too small, a sharp spike in the input data may kill the significant peak.

my $bam_file_vector = $arguments_tab[0];
my $input_file_vector = $arguments_tab[1];
my $experiment_name = $arguments_tab[2];
my $output_folder = $arguments_tab[3];
my $format = $arguments_tab[4];
my $genome_size = $arguments_tab[5];
my $tag_size = $arguments_tab[6];
my $qvalue = $arguments_tab[7];
my $elongation_size = $arguments_tab[8];
my $artefact_threshold = $arguments_tab[9];

if(! -d $output_folder)
{
	mkdir($output_folder) or die "Could not create $output_folder: $!";
}

my $output_folder_nomodel_broad = join('', $output_folder, "no_model_broad/");
my $output_folder_nomodel = join('', $output_folder, "no_model/");
my $output_folder_modelBased = join('', $output_folder, "model_based/");

if(! -d $output_folder_nomodel_broad)
{
	mkdir($output_folder_nomodel_broad) or die "Could not create $output_folder_nomodel_broad: $!";
}

if(! -d $output_folder_nomodel)
{
	mkdir($output_folder_nomodel) or die "Could not create $output_folder_nomodel: $!";
}

if(! -d $output_folder_modelBased)
{
	mkdir($output_folder_modelBased) or die "Could not create $output_folder_modelBased: $!";
}


print "This is job number $file_lineNumber\n";

if($input_file_vector ne "NA")
{
	print "---- Creating nomodel and broad\n";
	print "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_nomodel_broad -f $format -g $genome_size -s $tag_size --nomodel --extsize $elongation_size --keep-dup $artefact_threshold --broad --broad-cutoff $qvalue\n\n";
	my $commandToLoad = "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_nomodel_broad -f $format -g $genome_size -s $tag_size --nomodel --extsize $elongation_size --keep-dup $artefact_threshold --broad --broad-cutoff $qvalue";
	system($commandToLoad);
	
	print "---- Creating nomodel wihtout broad\n";
	print "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_nomodel -f $format -g $genome_size -s $tag_size -q $qvalue --nomodel --extsize $elongation_size --keep-dup $artefact_threshold\n\n";
	$commandToLoad = "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_nomodel -f $format -g $genome_size -s $tag_size -q $qvalue --nomodel --extsize $elongation_size --keep-dup $artefact_threshold";
	system($commandToLoad);
	
	print "---- Creating model based\n";
	print "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_modelBased -f $format -g $genome_size -s $tag_size -q $qvalue --keep-dup $artefact_threshold\n\n";
	$commandToLoad = "macs2 callpeak -t $bam_file_vector -c $input_file_vector -n $experiment_name --outdir $output_folder_modelBased -f $format -g $genome_size -s $tag_size -q $qvalue --keep-dup $artefact_threshold";
	system($commandToLoad);	
}else{
	print "---- Creating nomodel and broad\n";
	print "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_nomodel_broad -f $format -g $genome_size -s $tag_size --nomodel --extsize $elongation_size --keep-dup $artefact_threshold --broad --broad-cutoff $qvalue\n\n";
	my $commandToLoad = "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_nomodel_broad -f $format -g $genome_size -s $tag_size --nomodel --extsize $elongation_size --keep-dup $artefact_threshold --broad --broad-cutoff $qvalue";
	system($commandToLoad);
	
	print "---- Creating nomodel wihtout broad\n";
	print "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_nomodel -f $format -g $genome_size -s $tag_size -q $qvalue --nomodel --extsize $elongation_size --keep-dup $artefact_threshold\n\n";
	$commandToLoad = "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_nomodel -f $format -g $genome_size -s $tag_size -q $qvalue --nomodel --extsize $elongation_size --keep-dup $artefact_threshold";
	system($commandToLoad);
	
	print "---- Creating model based\n";
	print "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_modelBased -f $format -g $genome_size -s $tag_size -q $qvalue --keep-dup $artefact_threshold\n\n";
	$commandToLoad = "macs2 callpeak -t $bam_file_vector -n $experiment_name --outdir $output_folder_modelBased -f $format -g $genome_size -s $tag_size -q $qvalue --keep-dup $artefact_threshold";
	system($commandToLoad);
	
	
}





