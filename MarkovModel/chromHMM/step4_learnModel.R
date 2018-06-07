###########
# After having computed the binarized bed files, this script learn the HMM and output the different states.
# Descostes Dec 2016
###########

library("NGSprofiling");
library("Rargs");



################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--outfileID"]] <- list(variableName="outfileID", numeric=F, mandatory=T, description="The string outfileID is included in the file names of all the output files.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="Single string indicating the version of the genome. ex: hg19, mm10");
paramsDefinition[["--nb_cpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Nb of CPU used for parallelization.");
paramsDefinition[["--inputdir"]] <- list(variableName="inputdir", numeric=F, mandatory=T, description="Path to the folder containing input binarized bed files that were defined in step3.");
paramsDefinition[["--outputdir"]] <- list(variableName="outputdir", numeric=F, mandatory=T, description="Path to the folder where result files will be written.");
paramsDefinition[["--numstates"]] <- list(variableName="numstates", numeric=T, mandatory=T, description="Number of states given to the model.");




#Optional argument
paramsDefinition[["--segmentationSize"]] <- list(variableName="segmentation_size", numeric=T, mandatory=F, description="Should be the same than the one used in step3. Size of the sliding window used to perform segmentation of the genome. Default: 200 bp", default=200);
paramsDefinition[["--convergedelta"]] <- list(variableName="convergedelta", numeric=T, mandatory=F, description="The threshold on the change on the estimated log likelihood that if it falls below this value, then parameter training will terminate. If this value is less than 0 then it is not used as part of the stopping criteria. The default value for this parameter is 0.001.", default=0.001);
paramsDefinition[["--informationsmooth"]] <- list(variableName="informationsmooth", numeric=T, mandatory=F, description="A smoothing constant away from 0 for all parameters in the information based initialization. This option is ignored if random or load are selected for the initialization method. The default value of this parameter is 0.02.", default=0.02);
paramsDefinition[["--type_initialization"]] <- list(variableName="type_initialization", numeric=F, mandatory=F, description="possible values are information|random|load, default is information", default="information");
paramsDefinition[["--maxiterations"]] <- list(variableName="maxiterations", numeric=T, mandatory=F, description="This option specifies the maximum number of iterations over all the input data in the training. By default this is set to 200.", default=200);
paramsDefinition[["--stateordering"]] <- list(variableName="stateordering", numeric=F, mandatory=F, description="possible values are emission|transition. default: emission", default="emission");



#outfileID <- "model1";  # The string outfileID is included in the file names of all the output files.
#genome_version <- "hg19";
#nb_cpu <- 1;
#inputdir <- "/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/chromHMM/modele1/binary_files/";
#outputdir <- "/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/chromHMM/modele1/learned_states/"
#numstates <- 5;
#
##options
#
#segmentation_size <- 200; #Size of the sliding window used to perform segmentation of the genome. Default: 200 bp
#convergedelta <- 0.001; #The threshold on the change on the estimated log likelihood that if it falls below this value, then parameter training will terminate. If this value is less than 0 then it is not used as part of the stopping criteria. The default value for this parameter is 0.001.
#informationsmooth <- 0.02; #A smoothing constant away from 0 for all parameters in the information based initialization. This option is ignored if random or load are selected for the initialization method. The default value of this parameter is 0.02.
#type_initialization <- "information"; # possible values are information|random|load, default is information
#maxiterations <- 200; #This option specifies the maximum number of iterations over all the input data in the training. By default this is set to 200.
#stateordering <- "emission" #possible values are emission|transition. default: emission

# [-e loadsmoothemission] only used if "init" argument is load
# -holdcolumnorder
#-m modelinitialfile -This specifies the model file containing the initial parameters which can then be used with the load option
#-t loadsmoothtransition -This parameter is only applicable if the load option is selected for the init parameter. This parameter controls the smoothing away from 0 when loading a model. 

################



##############
# MAIN
##############



# Retreives the parameters
getParams(paramsDefinition);

checkingOutputFolder(outputdir);


valid_genomes <- c("hg18", "hg19", "hg38", "mm9", "mm10", "rn5", "rn6", "danRer7", "danRer10", "dm3", "dm6", "ce6", "ce10"); 

if(is.na(match(genome_version, valid_genomes)))
{
	stop("\n The name of the reference genome that you entered '", genome_version, "' is either no correct or not covered by chromHMM\n\n");
}


if(type_initialization == "load")
{
	stop("\n Script has to be modified to use an init value 'load', add the '-e' option to the script and create a if condition\n\n");
}

chromosomelengthfile <- paste("$CHROMHMM_ROOT/CHROMSIZES/", genome_version, ".txt", sep="");

		
command_to_load <- paste("java -jar $CHROMHMM_ROOT/ChromHMM.jar LearnModel -b", segmentation_size, "-d", convergedelta, "-h", informationsmooth, "-i", outfileID, "-init", type_initialization, "-l", chromosomelengthfile, "-p", nb_cpu, "-printposterior -printstatebyline -r", maxiterations, "-stateordering", stateordering, inputdir, outputdir, numstates, genome_version, sep=" ");


cat("Loading the command: ", command_to_load, "\n");

system(command_to_load);