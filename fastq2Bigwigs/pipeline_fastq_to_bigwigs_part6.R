############
# This script organizes the result of pasha and load the input subtraction (if chipseq = TRUE).
# Descostes August 2017
############




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
paramsDefinition[["--spikein"]]=list(variableName="spikein", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiments are spiked in.", postConversion=as.logical);
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="Name of the analysis. This will be used to create jobs names on the cluster.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
#chipseq <- TRUE;
#spikein <- TRUE;
#species <- "human";
#genome <- "hg19";

################



##############
# MAIN
##############

# Retreives the parameters
getParams(paramsDefinition);

if(single_end){
    wig_vec <- list.files(paste0(fastq_files_folder, "/bam_files/wig_files_SE/AllReads"), "wig", full.names = TRUE)
}else{
    wig_vec <- list.files(paste0(fastq_files_folder, "/bam_files/wig_files_PE/AllReads"), "wig", full.names = TRUE)
}

if(isTRUE(all.equal(length(wig_vec), 0)))
    stop("Wig files not found, problem in the pasha processing.\n");



## Organize files and load input subtraction if data are chip-seq experiments
setwd(paste0(fastq_files_folder, "/tmp"));
system("mkdir pasha_log");
Sys.sleep(3);
if(single_end) setwd(paste0(fastq_files_folder, "/bam_files/wig_files_SE")) else setwd(paste0(fastq_files_folder, "/bam_files/wig_files_PE")) 
system(paste0("mv ", fastq_files_folder, if(single_end) "/bam_files/wig_files_SE/*log" else "/bam_files/wig_files_PE/*log", " ", fastq_files_folder, "/tmp/pasha_log"));
system(paste0("mv ", fastq_files_folder, if(single_end) "/bam_files/wig_files_SE/AllReads/log_report" else "/bam_files/wig_files_PE/AllReads/log_report", " ", fastq_files_folder, "/tmp/pasha_log"));

setwd(paste0(fastq_files_folder));
system("mkdir wig_files");
Sys.sleep(3);
system(paste0("mv ", fastq_files_folder, if(single_end) "/bam_files/wig_files_SE/AllReads/*wig" else "/bam_files/wig_files_PE/AllReads/*wig", " ", fastq_files_folder, "/wig_files"));
system(paste0("rm -r ", fastq_files_folder, if(single_end) "/bam_files/wig_files_SE/" else "/bam_files/wig_files_PE/"));
Sys.sleep(3);



if(chipseq && !spikein) ## If chip-seq data and no spike-in, the input subtraction is performed here
{
    wig_vec <- list.files(paste0(fastq_files_folder, "/wig_files"), "wig", full.names = TRUE);
    input_file <- wig_vec[grep("nput", wig_vec)];
    wig_vec <- wig_vec[-which(wig_vec == input_file)];
    
    setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
    system("mkdir scale_and_subtract_input");
    Sys.sleep(3);
    setwd(paste0(fastq_files_folder, "/tmp/scripts/scale_and_subtract_input"));
    
    write(paste(input_file, 50, "NULL", "FALSE", "FALSE", "TRUE", sep=";"), file = "scaleInput.conf", ncolumns=1);
    
    write(paste(wig_vec, 50, input_file, "FALSE", "TRUE", "TRUE", sep=";"), file = "subtractInput.conf", ncolumns=1);
    line_end <- length(wig_vec);

    system(paste0("/ifs/home/descon01/cluster/scripts/executable/scale_and_subtract scaleInput_", analysis_name, " 1 1 scaleInput.conf"));
    
    command_sub <- paste0("/ifs/home/descon01/cluster/scripts/executable/scale_and_subtract_withHold scaleInput_", analysis_name, " inputSub_", analysis_name, " 1 ", line_end, " subtractInput.conf");
    system(command_sub);
}else if(chipseq && spikein){  ## The first step of the spike in normalization is to create indexes for bam files
    
    bam_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), "_SORTED_PICARD_COOR.bam", full.names = TRUE);
    
    setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
    system("mkdir spikeIn_process");
    Sys.sleep(3);
    setwd(paste0(fastq_files_folder, "/tmp/scripts/spikeIn_process"));
    
    write(bam_files_vec, file="createBamIndex.conf", ncolumns=1);
    line_end <- length(bam_files_vec);
    
    command_bamIndex <- paste0("/ifs/home/descon01/cluster/scripts/executable/create_bam_index BamIndex_", analysis_name, " 1 ", line_end, " createBamIndex.conf");
    system(command_bamIndex);
    
}  



## ------------------------ Launch the next step

setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_6 <- "#!/bin/bash";
command_part_6 <- c(command_part_6, "#$ -cwd");
command_part_6 <- c(command_part_6, "module load r/3.2.0");
command_part_6 <- c(command_part_6, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part7.R --fastqFilesFolder ", fastq_files_folder, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --spikein ", spikein, " --analysisName ", analysis_name));

write(command_part_6, file="load_part7.sh", ncolumns=1);

if(chipseq && !spikein)
{
    command_part6 <- paste("qsub -q all.q -hold_jid inputSub_", analysis_name, " load_part7.sh", sep="");
}else if(chipseq && spikein){
    command_part6 <- paste("qsub -q all.q -hold_jid BamIndex_", analysis_name, " load_part7.sh", sep="");
}else{
    command_part6 <- paste("qsub -q all.q load_part7.sh", sep="");
}

system(command_part6);






