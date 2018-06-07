###############
# This script performs peak calling with GEM (GPS) software.
# Descostes march 2017
###############

#with java/1.7

library("NGSprofiling");
library("Rargs");


################
# PARAMETERS
################


#parameters defined from the command line using RIO
paramsDefinition <- list();

paramsDefinition[["--experimentBed"]] <- list(variableName="experiment_bed", numeric=F, mandatory=T, description="Single file path to the experiment bed file.");
paramsDefinition[["--inputBed"]] <- list(variableName="input_bed", numeric=F, mandatory=T, description="Single file path to the input bed file.");
paramsDefinition[["--genomeChromsizeFile"]] <- list(variableName="genome_chromsize_file", numeric=F, mandatory=T, description="File provided with the package and containing chrom sizes ot the corresponding genome.");
paramsDefinition[["--qValue"]] <- list(variableName="q_value", numeric=T, mandatory=T, description="Q-value threshold, will be converted to -log10 by the script.");
paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Number of cpu to use for normalization.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Single string giving the path to the output folder.");
paramsDefinition[["--expname"]] <- list(variableName="expname", numeric=F, mandatory=T, description="Single string giving the name of the experiment.");


#experiment_bed <- "/ifs/home/descon01/data/data_december2016/bed_from_bam/james_K27M_chipseq_spikedin/EZH2_H33WT_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bed";
#input_bed <- "/ifs/home/descon01/data/data_december2016/bed_from_bam/james_K27M_chipseq_spikedin/input_H33WT_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bed";
#genome_chromsize_file <- "/ifs/home/descon01/programs/gem/hg19.chrom.sizes";
#q_value <- 0.01;
#nb_cpu <- 1;
#output_folder <- "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_GEM/";
#expname <- "EZH2_H33WT";

################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


checkingOutputFolder(output_folder);


if(length(experiment_bed) > 1)
{
	stop("\n This script treats only one experiment at once\n");
}

minus_log10_qvalue <- -log10(q_value);
output_file_name <- paste(output_folder, expname, sep="");

#GEM can be activated by giving a genome sequence (--genome) and using any one of the following command line options:
		
#--k: the length of the k-mers
#--k_min and --k_max: the range for the length of k-mers
#--seed: the seed k-mer to jump start k-mer set motif discovery. The length of the seed k-mer will be used to set k. 
#If these three options are not used, GEM will just run GPS and stop.


command <- paste("java -Xmx15G -jar /ifs/home/descon01/programs/gem/gem.jar --d /ifs/home/descon01/programs/gem/Read_Distribution_default.txt --exptX ", experiment_bed, " --ctrlX ", input_bed, " --g ", genome_chromsize_file, " --q ", minus_log10_qvalue, " --t ", nb_cpu, " --out ", output_file_name, sep="");

cat("Loading ", command, "\n");
system(command);
