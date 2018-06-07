#############
## This script moves the bam files to a bam folder for chipseq and then performs sorting.
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
#single_end <- TRUE;
#fragment_size <- 150;
#wig_fixed <- TRUE;
#wig_variable <- FALSE;
#spikein <- TRUE;
#species <- "human";
#genome <- "hg19";

################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

Sys.sleep(300);


if(chipseq)
{
    setwd(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"));
    
    if(!isTRUE(all.equal(species, "ant"))){
        bam_files_vec <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), ".bam", full.names = TRUE);
        log_files_vec <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), "_LOG", full.names = TRUE);
    }else{
        bam_files_vec <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), ".bam", full.names = TRUE, recursive=TRUE)
        log_files_vec <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), "final.out", full.names = TRUE, recursive = TRUE)
    }
    
    
    if(isTRUE(all.equal(length(bam_files_vec), 0)))
        stop("Bam files not found, problem in the conversion of sam files.\n");
    
    bam_folder <- paste0(fastq_files_folder, "/bam_files");
    system(paste0("mkdir ", bam_folder));
    
    Sys.sleep(1);
    sapply(paste0("mv ", bam_files_vec, " ", fastq_files_folder, "/bam_files"), function(x){system(x)});
    sapply(paste0("mv ", log_files_vec, " ", fastq_files_folder, "/bam_files"), function(x){system(x)});
}


## Performing bam files sorting with PICARD, if STAR was used, the bam files are already sorted

if(!isTRUE(all.equal(species, "ant"))){
    bam_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), ".bam", full.names = TRUE);
    
    if(isTRUE(all.equal(length(bam_files_vec), 0)))
        stop("Bam files not found, problem in the conversion of sam files or in alignment of rna-seq data.\n");
    
    setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
    system("mkdir sort_bam");
    Sys.sleep(1);
    setwd(paste0(fastq_files_folder, "/tmp/scripts/sort_bam"));
    write(bam_files_vec, file="sort_bam.conf", ncolumns=1);
    Sys.sleep(5);
    line_end <- length(bam_files_vec);
    
    command_sort <- paste0("/ifs/home/descon01/cluster/scripts/executable/sort_bam_bycoor sort_bam_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " sort_bam.conf");
    
    system(command_sort);
}


## ------------------------ Launch the next step


setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_4 <- "#!/bin/bash";
command_part_4 <- c(command_part_4, "#$ -cwd");
command_part_4 <- c(command_part_4, "module load r/3.2.0");
command_part_4 <- c(command_part_4, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part5.R --fastqFilesFolder ", fastq_files_folder, " --chipseq ", chipseq, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --fragmentSize ", fragment_size, " --wigFixed ", wig_fixed, " --wigVariable ", wig_variable, " --spikein ", spikein, " --analysisName ", analysis_name));

write(command_part_4, file="load_part5.sh", ncolumns=1);

command_part4 <- paste("qsub -q all.q -hold_jid sort_bam_", analysis_name, " load_part5.sh", sep="");

system(command_part4);
