#############
## This script performs the first steps of data processing which are quality control and quality filtering.
## Descostes August 2017
#############


library(Rargs);


################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


paramsDefinition[["--fastqFilesFolder"]] <- list(variableName="fastq_files_folder", numeric=F, mandatory=T, description="Folder paths to the fastq files.");
paramsDefinition[["--chipseq"]]=list(variableName="chipseq", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiment is of type chipseq", postConversion=as.logical);
paramsDefinition[["--genome"]] <- list(variableName="genome", numeric=F, mandatory=T, description="Genome version to align the data on. Currently hg19 and mm10 are supported.");
paramsDefinition[["--species"]] <- list(variableName="species", numeric=F, mandatory=T, description="Name of species. Currently human and mouse are supported.");
paramsDefinition[["--singleEnd"]]=list(variableName="single_end", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiment is single or paired ended. Only single end is currently supported.", postConversion=as.logical);
paramsDefinition[["--fragmentSize"]] <- list(variableName="fragment_size", numeric=T, mandatory=T, description="Size of the sonicated fragments. Should be 0 for rna-seq.");
paramsDefinition[["--wigFixed"]]=list(variableName="wig_fixed", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if wiggle files in fixed steps should be output.", postConversion=as.logical);
paramsDefinition[["--wigVariable"]]=list(variableName="wig_variable", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if wiggle files in variable steps should be output.", postConversion=as.logical);
paramsDefinition[["--spikein"]]=list(variableName="spikein", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiments are spiked in.", postConversion=as.logical);
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="Name of the analysis. This will be used to create jobs names on the cluster.");
paramsDefinition[["--minCpuNumber"]] <- list(variableName="min_cpu_number", numeric=T, mandatory=T, description="For processes using parallelization, minimum number of cpu that should be available to run the script.");
paramsDefinition[["--maxCpuNumber"]] <- list(variableName="max_cpu_number", numeric=T, mandatory=T, description="For processes using parallelization, maximum number of cpu that could be used to run the script.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
#chipseq <- TRUE;
#genome <- "hg19";
#species <- "human";
#single_end <- TRUE;
#fragment_size <- 150;
#wig_fixed <- TRUE;
#wig_variable <- FALSE;
#spikein <- TRUE;

#---------------------------------------- RNA-seq


#fastq_files_folder <- "/ifs/home/descon01/analysis/core_facility_lab/rna-seq";
#chipseq <- FALSE;
#genome <- "mm10";
#species <- "mouse";
#single_end <- TRUE;
#fragment_size <- 150;
#wig_fixed <- TRUE;
#wig_variable <- FALSE;
#spikein <- FALSE;


fastq_files_folder <- "/ifs/home/descon01/data/data_march2018/fasteq/comzit_workerfatbody"
chipseq <- TRUE
genome <- "hsal8_5"
species <- "ant"
single_end <- FALSE
fragment_size <- 150
wig_fixed <- TRUE
wig_variable <- FALSE
spikein <- FALSE
analysis_name <- "comzit_workerfatbody"
min_cpu_number <- 10
max_cpu_number <- 20



################



##############
# MAIN
##############

# Retreives the parameters
getParams(paramsDefinition);


if(!dir.exists(fastq_files_folder)){
    stop(fastq_files_folder, " does not exist or is not a valid folder\n");
}

fastq_files <- list.files(fastq_files_folder, ".fastq", full.names = TRUE);

if(isTRUE(all.equal(length(fastq_files), 0))){
    stop("No fastq files were found in the folder ", fastq_files_folder, "\n");
}

if(!isTRUE(all.equal(length(grep("gz", fastq_files)),0))){
    stop("Your files shound not be compressed, please uncompress before processing\n");
}

if(!(isTRUE(all.equal(species, "human")) || isTRUE(all.equal(species, "mouse")) || isTRUE(all.equal(species, "ant")))){
    stop("The species should be ant, human or mouse, other organims supported on demand.\n");
}

if(chipseq && isTRUE(all.equal(length(grep("nput", fastq_files)), 0)))
{
    stop("The input file was not found in your folder, the file name should contain 'nput' characters\n");
}


if(chipseq && single_end && !isTRUE(all.equal(length(grep("nput", fastq_files)), 1)))
{
    stop("Only one input file should be provided. If your fastq files represent experiments in several conditions, create as many folders as conditions and process files separately.\n");
}

if(chipseq && !single_end && !isTRUE(all.equal(length(grep("nput", fastq_files)), 2)))
{
    stop("Only one input with each pair file should be provided. You should therefore provide two files for the input. If your fastq files represent experiments in several conditions, create as many folders as conditions and process files separately.\n");
}

if(!chipseq && spikein){
    stop("The spikein rna-seq is not handled by this pipeline\n");
}

if(!chipseq && !single_end)
    stop("Paired-end rna-seq is not yet supported\n")

if(!single_end && !isTRUE(all.equal(length(fastq_files)%%2, 0)))
{
    stop("If your experiments are paired-ended, the number of files should be even\n");
}

if(!wig_fixed && !wig_variable)
{
    stop("At least one type of wiggle file should be output, fixed or variable.")
}

if(!single_end && !isTRUE(all.equal(length(grep("_R.+\\.fastq", fastq_files)), length(fastq_files))))
    stop("Each experiment files should end by _R1.fastq or _R2.fastq to indicate each pair.")

if(chipseq && spikein && isTRUE(all.equal(species, "ant"))){
    stop("Spike-in is not supported for ant")
}

if(!chipseq && isTRUE(all.equal(species, "ant"))){
    stop("RNA-seq is not yet supported for ant")
}

## ------------- Launching quality control

## Creating folders to launch scripts
setwd(fastq_files_folder);
system("mkdir tmp");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp"));
system("mkdir scripts");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp/scripts"));
system("mkdir quality_control");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp/scripts/quality_control"));

## Preparing files to launch quality control
write(fastq_files, file="quality_control.conf", ncolumns=1);
line_end <- length(fastq_files)

command_qc <- paste0("/ifs/home/descon01/cluster/scripts/executable/qc_launcher quality_control_", analysis_name, " 1 ", line_end, " quality_control.conf");
system(command_qc);




##--------------- Launching quality filtering


setwd(paste0(fastq_files_folder, "/tmp/scripts"));
system("mkdir quality_filtering");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp/scripts/quality_filtering"));

# Preparing file to launch quality filtering

if(single_end){
    write(paste(fastq_files, 80, 25, sep=";"), file="quality_filtering.conf", ncolumns=1);
    command_qf <- paste0("/ifs/home/descon01/cluster/scripts/executable/qc_filtering_singleEnd quality_control_", analysis_name, " quality_filtering_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " quality_filtering.conf");
}else{
    
    # Organizing fastq files by pairs
    fastq_factors <-as.factor(unlist(lapply(strsplit(fastq_files, "_R"), "[[", 1)))
    fastq_grouped_list <- split(fastq_files, fastq_factors)
    to_write_qfilPaired <- as.character(unlist(lapply(fastq_grouped_list, function(x){return(paste(x[1], x[2], 80, 25, sep=";"))})))
    write(to_write_qfilPaired, file="quality_filtering.conf", ncolumns=1);
    Sys.sleep(10)
    line_end <- line_end/2
    command_qf <- paste0("/ifs/home/descon01/cluster/scripts/executable/qc_filtering_pairedEnd quality_control_", analysis_name, " quality_filtering_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " quality_filtering.conf");
}

system(command_qf);

##------------------ Launching the next step of pipeline

setwd(paste0(fastq_files_folder, "/tmp/scripts"));
system("mkdir steps_transition");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));


command_part_1 <- "#!/bin/bash";
command_part_1 <- c(command_part_1, "#$ -cwd");
command_part_1 <- c(command_part_1, "module load r/3.2.0");
command_part_1 <- c(command_part_1, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part2.R --fastqFilesFolder ", fastq_files_folder, " --chipseq ", chipseq, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --fragmentSize ", fragment_size, " --wigFixed ", wig_fixed, " --wigVariable ", wig_variable, " --spikein ", spikein, " --analysisName ", analysis_name, " --minCpuNumber ", min_cpu_number, " --maxCpuNumber ", max_cpu_number));

write(command_part_1, file="load_part2.sh", ncolumns=1);
command_part1 <- paste("qsub -q all.q -hold_jid quality_filtering_", analysis_name, " load_part2.sh", sep="");
system(command_part1);



