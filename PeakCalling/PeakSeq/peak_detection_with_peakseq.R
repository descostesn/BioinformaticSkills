############
# This script first creates a configuration file and then run peakseq peak detection.
# Descostes March 2017
############

library("Rargs");
library("NGSprofiling");

################
# PARAMETERS
################


#parameters defined from the command line using RIO
paramsDefinition <- list();

paramsDefinition[["--expId"]] <- list(variableName="exp_id", numeric=F, mandatory=T, description="Experiment id is used as a prefix to the output file name..");
paramsDefinition[["--tagLength"]] <- list(variableName="tag_length", numeric=T, mandatory=T, description="Average read size.");
paramsDefinition[["--targetFdr"]] <- list(variableName="target_fdr", numeric=T, mandatory=T, description="Target FDR in the simulations.");
paramsDefinition[["--numberSimulation"]] <- list(variableName="number_simulation", numeric=T, mandatory=T, description="Number of simulation.");
paramsDefinition[["--minInterpeakDistance"]] <- list(variableName="min_interpeak_distance", numeric=T, mandatory=T, description="Distance minimum between two detected peaks.");
paramsDefinition[["--mappabilityMapFile"]] <- list(variableName="mappability_map_file", numeric=F, mandatory=T, description="Mappability file that includes the uniquely mappable number of nucleotides per window for each chromosome.");
paramsDefinition[["--chipseqReadsFolder"]] <- list(variableName="chipseq_reads_folder", numeric=F, mandatory=T, description="Path to the folder containing the pre-processed chip-seq data of the considered experiment.");
paramsDefinition[["--inputReadsFolder"]] <- list(variableName="input_reads_folder", numeric=F, mandatory=T, description="Path to the folder containing the pre-processed chip-seq data of the control experiment.");
paramsDefinition[["--Qvalue"]] <- list(variableName="q_value", numeric=T, mandatory=T, description="Q-value threshold applied on the final set of peaks.");
paramsDefinition[["--backgroundModelMode"]] <- list(variableName="background_model_mode", numeric=F, mandatory=T, description="Can be 'Simulated' or 'Poisson'.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the output folder.");


################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

checkingOutputFolder(output_folder);

config_file <- paste("# Experiment id is used as a prefix to the output file name.\n",
	  "Experiment_id ", exp_id, "\n",
	  "# Enrichment fragment length For tag extension, this is the value of average fragment length.\n",
	  "Enrichment_mapped_fragment_length ", tag_length, "\n",
	  "# Target FDR in the simulations.\n",
	  "target_FDR ", target_fdr, "\n",
	  "# Number of simulations performed while estimating the putative peaks.\n",
	  "N_Simulations ", number_simulation, "\n",
	  "# Minimum distance between consecutive peaks\n",
      "Minimum_interpeak_distance ", min_interpeak_distance, "\n",
      "# Mappability file that includes the uniquely mappable number of nucleotides per window for each chromosome.\n",
      "Mappability_map_file ", mappability_map_file, "\n",
      "# The directory that contains the preprocessed ChIP-Seq reads, can specify multiple directories to pool reads from multiple source (e.g. replicates)\n",
      "ChIP_Seq_reads_data_dirs ", chipseq_reads_folder, "\n",
      "# The directory that contains the preprocessed Input (control) experiment reads. (Multiple directories allowed)\n",
      "Input_reads_data_dirs ", input_reads_folder, "\n",
	  "# Seed for pseudo-random number generator. This is necessary for simulated background option (specified below).\n",
      "#Simulation_seed 1234567\n",
	  "# Q-value threshold applied on the final set of peaks.\n",
	  "max_Qvalue ", q_value, "\n",
      "# There are currently two models for simulating the background for threshold selection\n
	   # Simulated background is the simulation based method that is explained in the PeakSeq paper.\n
       # Poisson background uses a simple Poisson background with mean estimated from the read statistics. This option is still experimental but it is much faster than the simulated background option.\n
       # Background_model Poisson\n",
      "Background_model ", background_model_mode, "\n", sep="");
		
write(config_file, file=paste(output_folder, "param.config", sep=""), ncolumns=1);


command_to_load <- paste("/ifs/home/descon01/programs/PeakSeq/bin/PeakSeq -peak_select ", output_folder, "param.config", sep="");
system(command_to_load);
		






