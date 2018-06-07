############
# If spikein = TRUE, this script performs all steps for normalization.
# Descostes August 2017
############


library(Rargs);


################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


paramsDefinition[["--fastqFilesFolder"]] <- list(variableName="fastq_files_folder", numeric=F, mandatory=T, description="Folder paths to the fastq files.");
paramsDefinition[["--genome"]] <- list(variableName="genome", numeric=F, mandatory=T, description="Genome version to align the data on. Currently hg19 and mm10 are supported.");
paramsDefinition[["--species"]] <- list(variableName="species", numeric=F, mandatory=T, description="Name of species. Currently human and mouse are supported.");
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="Name of the analysis. This will be used to create jobs names on the cluster.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
#species <- "human";
#genome <- "hg19";



#fastq_files_folder <- "/ifs/home/descon01/analysis/core_facility_lab/chipseq/spikedin/mouse/";
#spikein <- TRUE;
#species <- "mouse";
#genome <- "mm10";
#analysis_name <- "chipseqspikedinmouse";




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);



##
## Perform RPM normalization
##

wig_vec <- list.files(paste0(fastq_files_folder, "/wig_files"), "wig", full.names = TRUE);

factors_file_endo <- paste0(fastq_files_folder, "/bam_files/scores_", genome, ".txt");

if(!file.exists(factors_file_endo)){
    stop("The endo scores file was not found, problem in scaling factors computing\n");
}

factors_endo <- readLines(factors_file_endo);
expnames_file_factors <- unlist(strsplit(unlist(lapply(strsplit(factors_endo,":"), "[", 1)), paste("_", genome, sep="")));
factors <- as.numeric(unlist(lapply(strsplit(factors_endo,":"), "[", 2)));

## Order wig vec
index_sort_wig <- as.numeric(sapply(expnames_file_factors, function(x,wig){return(grep(x,wig))}, wig_vec));
wig_vec <- wig_vec[index_sort_wig];

setwd(paste0(fastq_files_folder, "/tmp/scripts/spikeIn_process"));

param_conf_RPM <- paste(wig_vec, factors, 50, sep=";");
write(param_conf_RPM, file= "applyRPM.conf",ncolumns=1);

line_end <- length(param_conf_RPM);

command_RPM <- paste0("/ifs/home/descon01/cluster/scripts/executable/apply_RPM RPM_", analysis_name, " 1 ", line_end, " applyRPM.conf");
system(command_RPM);

##
## Perform input subtraction
##


wig_vec <- paste0(unlist(strsplit(wig_vec, "\\.wig")), "-RPM.wig");
input_file <- wig_vec[grep("nput", wig_vec)];
wig_vec <- wig_vec[-which(wig_vec == input_file)];

write(paste(wig_vec, 50, input_file, "FALSE", "TRUE", "FALSE", sep=";"), file = "subtractInput.conf", ncolumns=1);
line_end <- length(wig_vec);
command_sub <- paste0("/ifs/home/descon01/cluster/scripts/executable/scale_and_subtract_withHold RPM_", analysis_name, " inputSubSpike_", analysis_name, " 1 ", line_end, " subtractInput.conf");
system(command_sub);

##
## Reverse RPM normalization
##

wig_vec <- paste(unlist(strsplit(wig_vec, "\\.wig")), "_BGSub.wig",sep="");
index_remove_input <- as.numeric(sapply(expnames_file_factors, function(x,input){return(grep(x,input))}, input_file));
index_remove_input <- index_remove_input[which(!is.na(index_remove_input))];
factors <- factors[-index_remove_input];

param_reverse <- paste(wig_vec, factors, 50, sep=";"); 
write(param_reverse, file= "reverse_RPM.conf", ncolumns=1);

command_reverse <- paste0("/ifs/home/descon01/cluster/scripts/executable/reverse_RPM_withhold inputSubSpike_", analysis_name, " reverse_", analysis_name, " 1 ", line_end, " reverse_RPM.conf");
system(command_reverse);



##
## Apply exo scaling
##

wig_vec <- paste(unlist(strsplit(wig_vec, "\\.wig")), "-scaleReverse.wig",sep="");
factors_file_exo <- paste0(fastq_files_folder, "/bam_files/scores_dm3.txt");

if(!file.exists(factors_file_exo)){
    stop("The exo scores file was not found, problem in scaling factors computing\n");
}

factors_exo <- readLines(factors_file_exo);
expnames_file_factors <- unlist(strsplit(unlist(lapply(strsplit(factors_exo,":"), "[", 1)), "_dm3"));
factors <- as.numeric(unlist(lapply(strsplit(factors_exo,":"), "[", 2)));
index_remove_input <- as.numeric(sapply(expnames_file_factors, function(x,input){return(grep(x,input))}, input_file));
#index_remove_input <- index_remove_input[which(!is.na(index_remove_input))];
#factors <- factors[-index_remove_input];
#expnames_file_factors <- expnames_file_factors[-index_remove_input];
#index_sort_wig <- as.numeric(sapply(expnames_file_factors, function(x,wig){return(grep(x,wig))}, wig_vec));
#wig_vec <- wig_vec[index_sort_wig];

index_remove_input <- which(!is.na(index_remove_input));
factors <- factors[-index_remove_input];
expnames_file_factors <- expnames_file_factors[-index_remove_input];
index_sort_wig <- as.numeric(sapply(expnames_file_factors, function(x,wig){return(grep(x,wig))}, wig_vec));
wig_vec <- wig_vec[index_sort_wig];

param_exoscaling <- paste(wig_vec, factors, 50, sep=";"); 
write(param_exoscaling, file="exo_scaling.conf", ncolumns=1);

command_exoScaling <- paste0("/ifs/home/descon01/cluster/scripts/executable/apply_exoScaling_withhold reverse_", analysis_name, " exoScaling_", analysis_name, " 1 ", line_end, " exo_scaling.conf");
system(command_exoScaling);



#------------- launch next script

setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_8 <- "#!/bin/bash";
command_part_8 <- c(command_part_8, "#$ -cwd");
command_part_8 <- c(command_part_8, "module load r/3.2.0");
command_part_8 <- c(command_part_8, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part9.R --fastqFilesFolder ", fastq_files_folder, " --species ", species));

write(command_part_8, file="load_part9.sh", ncolumns=1);

command_part8 <- paste("qsub -q all.q -hold_jid exoScaling_", analysis_name, " load_part9.sh", sep="");
system(command_part8);







