#############
## This script converts the sam to bam for chipseq and move the bam files to a bam folder for rna-seq.
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
paramsDefinition[["--filteredFastqFiles"]] <- list(variableName="filtered_fastq_files", numeric=F, mandatory=T, description="Obtained from the previous step.");
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="Name of the analysis. This will be used to create jobs names on the cluster.");
paramsDefinition[["--minCpuNumber"]] <- list(variableName="min_cpu_number", numeric=T, mandatory=T, description="For processes using parallelization, minimum number of cpu that should be available to run the script.");
paramsDefinition[["--maxCpuNumber"]] <- list(variableName="max_cpu_number", numeric=T, mandatory=T, description="For processes using parallelization, maximum number of cpu that could be used to run the script.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
#chipseq <- TRUE;
#filtered_fastq_files <- OBTAINED FROM PREVIOUS SCRIPT
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


setwd(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"));

if(!isTRUE(all.equal(species, "ant"))){
    if(chipseq)
    {
        sam_files_vec <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), ".sam", full.names=TRUE);
        
        if(isTRUE(all.equal(length(sam_files_vec),0)))
            stop("The sam files were not found, there is a problem with bowtie\n");
        
        setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
        system("mkdir sam_to_bam");
        Sys.sleep(5);
        setwd(paste0(fastq_files_folder, "/tmp/scripts/sam_to_bam"));
        
        write(sam_files_vec, file="samBam.conf", ncolumns=1);
        
        line_end <- length(sam_files_vec);
        command_convert <- paste0("/ifs/home/descon01/cluster/scripts/executable/sam_to_bam sambam_", analysis_name, " 1 ", line_end, " samBam.conf");
        system(command_convert);
    }else{
        
        exp_name_vec <- unlist(lapply(strsplit(basename(filtered_fastq_files), "_filtered.fastq"),"[",1));
        setwd("./tophat_singleEnd_UNI");
        
        #Converting names
        sapply(exp_name_vec, function(x){system(paste0("mv ", x, "/accepted_hits.bam ", x, "/", x, ".bam"))});
        sapply(exp_name_vec, function(x){system(paste0("mv ", x, "/align_summary.txt ", x, "/", x, "_LOG.txt"))});
        
        bam_folder <- paste0(fastq_files_folder, "/bam_files");
        system(paste0("mkdir ", bam_folder));
        
        # Moving files to bam folder
        Sys.sleep(5);
        sapply(exp_name_vec, function(x){system(paste0("mv ", x, "/", x, ".bam ", bam_folder))});
        sapply(exp_name_vec, function(x){system(paste0("mv ", x, "/",x, "_LOG.txt ", bam_folder))});
    }
}


## ----------------------- Launching the next step


setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_3 <- "#!/bin/bash";
command_part_3 <- c(command_part_3, "#$ -cwd");
command_part_3 <- c(command_part_3, "module load r/3.2.0");
command_part_3 <- c(command_part_3, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part4.R --fastqFilesFolder ", fastq_files_folder, " --chipseq ", chipseq, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --fragmentSize ", fragment_size, " --wigFixed ", wig_fixed, " --wigVariable ", wig_variable, " --spikein ", spikein, " --analysisName ", analysis_name, " --minCpuNumber ", min_cpu_number, " --maxCpuNumber ", max_cpu_number));

write(command_part_3, file="load_part4.sh", ncolumns=1);

if(chipseq){
    command_part3 <- paste("qsub -q all.q -hold_jid sambam_", analysis_name, " load_part4.sh", sep="");
}else{
    command_part3 <- paste("qsub -q all.q load_part4.sh", sep="");
}

system(command_part3);



