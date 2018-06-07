#######
# This script binarizes bed files that were obtained from bam files in step1.
# Descostes December 2016
#######

library("NGSprofiling");
library("Rargs");



################
# PARAMETERS
################



#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--controlDirectory"]] <- list(variableName="control_directory", numeric=F, mandatory=T, description="Path to the folder containing input bed files that were defined in step2.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="Single string indicating the version of the genome. ex: hg19, mm10");
paramsDefinition[["--inputbeddir"]] <- list(variableName="inputbeddir", numeric=F, mandatory=T, description="Single string indicating the path to the bed files.");
paramsDefinition[["--cellmarkfiletable"]] <- list(variableName="cellmarkfiletable", numeric=F, mandatory=T, description="Path to the file created in step2.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the folder where binarized files will be written.");


#Optional argument
paramsDefinition[["--segmentationSize"]] <- list(variableName="segmentation_size", numeric=T, mandatory=F, description="Size of the sliding window used to perform segmentation of the genome. Default: 200 bp", default=200);
paramsDefinition[["--foldthresh"]] <- list(variableName="foldthresh", numeric=T, mandatory=F, description="This indicates a threshold for the fold enrichment over expected that must be met or exceeded by the observed count in a bin for a present call. By default this parameter value is 0 meaning effectively it is not used.", default=0);
paramsDefinition[["--signalthresh"]] <- list(variableName="signalthresh", numeric=T, mandatory=F, description="This indicates a threshold for the signal that must be met or exceeded by the observed count in a bin for a present call. This parameter can be useful when desiring to directly place a threshold on the signal. By default this parameter value is 0 meaning effectively it is not used.", default=0);
paramsDefinition[["--shift"]] <- list(variableName="shift", numeric=T, mandatory=F, description="Correspond to the elongation size determined by pasha. Take a median of all the values. default: 100.", default=100);
paramsDefinition[["--poissonthresh"]] <- list(variableName="poissonthresh", numeric=T, mandatory=F, description="This option specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter is 0.0001.", default=0.0001);
paramsDefinition[["--pseudocountcontrol"]] <- list(variableName="pseudocountcontrol", numeric=T, mandatory=F, description="An integer pseudocount that is uniformly added to every bin in the control data in order to smooth the control data from 0. The default value is 1.", default=1);
paramsDefinition[["--flankwidthcontrol"]] <- list(variableName="flankwidthcontrol", numeric=T, mandatory=F, description="This determines the extent of the spatial smoothing in computing the local enrichment for control reads. The local enrichment for control signal in the xth bin on the chromosome after adding pseudocountcontrol is computed based on the average control counts for all bins within x-w and x+w. If no controldir is specified, then this option is ignored. The default value is 5.", default=5);

#control_directory <- "/ifs/home/descon01/data/gary_paper_data/bed_files_from_bam/"; #Path to the folder containing input bed files that were defined in step2.
#genome_version <- "hg19";
#inputbeddir <- "/ifs/home/descon01/data/gary_paper_data/bed_files_from_bam/";
#cellmarkfiletable <- "/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/chromHMM/modele1/marks_ledgf_model1.txt";
#output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/chromHMM/modele1/binary_files/"; 
#
##optional
#
#segmentation_size <- 200; #Size of the sliding window used to perform segmentation of the genome. Default: 200 bp
#foldthresh <- 0;          #This indicates a threshold for the fold enrichment over expected that must be met or exceeded by the observed count in a bin for a present call. By default this parameter value is 0 meaning effectively it is not used.
#signalthresh <- 0; #This indicates a threshold for the signal that must be met or exceeded by the observed count in a bin for a present call. This parameter can be useful when desiring to directly place a threshold on the signal. By default this parameter value is 0 meaning effectively it is not used.
#shift <- 100;    # Correspond to the elongation size determined by pasha. Take a median of all the values. default: 100
#poissonthresh <- 0.0001; # This option specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter is 0.0001
#pseudocountcontrol <- 1; #An integer pseudocount that is uniformly added to every bin in the control data in order to smooth the control data from 0. The default value is 1.
#flankwidthcontrol <- 5; # This determines the extent of the spatial smoothing in computing the local enrichment for control reads. The local enrichment for control signal in the xth bin on the chromosome after adding pseudocountcontrol is computed based on the average control counts for all bins within x-w and x+w. If no controldir is specified, then this option is ignored. The default value is 5. 

################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


checkingOutputFolder(output_folder);


valid_genomes <- c("hg18", "hg19", "hg38", "mm9", "mm10", "rn5", "rn6", "danRer7", "danRer10", "dm3", "dm6", "ce6", "ce10"); 

if(is.na(match(genome_version, valid_genomes)))
{
	stop("\n The name of the reference genome that you entered '", genome_version, "' is either no correct or not covered by chromHMM\n\n");
}

outputcontroldir <- output_folder;
outputsignaldir <- output_folder;
chromosomelengthfile <- paste("$CHROMHMM_ROOT/CHROMSIZES/", genome_version, ".txt", sep="");



command_to_load <- paste("java -jar $CHROMHMM_ROOT/ChromHMM.jar BinarizeBed", "-b", segmentation_size, "-c", control_directory, "-f", foldthresh, "-g", signalthresh, "-n", shift, "-o", outputcontroldir,
"-p", poissonthresh, "-t", outputsignaldir, "-u", pseudocountcontrol, "-w", flankwidthcontrol, chromosomelengthfile, inputbeddir, cellmarkfiletable, output_folder, sep=" ");

cat("Loading the command: ", command_to_load, "\n");

system(command_to_load);