#############
## This script organizes the output files of the quality control and filtering. It then performs the alignment of the data.
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


fastq_files_folder <- "/ifs/home/descon01/data/data_march2018/fasteq/thelma"
chipseq <- TRUE
genome <- "mm10"
species <- "mouse"
single_end <- FALSE
fragment_size <- 150
wig_fixed <- TRUE
wig_variable <- FALSE
spikein <- TRUE
analysis_name <- "thelma"
min_cpu_number <- 1
max_cpu_number <- 20
#
#
#fastq_files_folder <- "/ifs/home/descon01/data/data_march2018/fasteq/comzit_gamergatefatbody"
#chipseq <- TRUE
#genome <- "hsal8_5"
#species <- "ant"
#single_end <- FALSE
#fragment_size <- 150
#wig_fixed <- TRUE
#wig_variable <- FALSE
#spikein <- FALSE
#analysis_name <- "comzit_gamergatefatbody"
#min_cpu_number <- 10
#max_cpu_number <- 20



################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

Sys.sleep(300);
options(scipen=999)

##------------------ Organize the output files of the quality control and filtering

setwd(fastq_files_folder);
system("mkdir quality_control_and_filtering");
Sys.sleep(1);
system("mv IlluQC_Filtered_files/ quality_control_and_filtering");
system("mv *html quality_control_and_filtering");
system("mv *zip quality_control_and_filtering");


##------------------ Performing alignment of the data

## Retrieving filtered fastq files

setwd(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files/"));
system("ls *fastq_filtered | sed -e \"p;s/\\.fastq_filtered/_filtered\\.fastq/\" | xargs -n2 mv");
Sys.sleep(1);
filtered_fastq_files <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), "_filtered.fastq", full.names = TRUE);
filtered_stat_files <- list.files(paste0(fastq_files_folder, "/quality_control_and_filtering/IlluQC_Filtered_files"), "fastq_stat", full.names = TRUE);

if(isTRUE(all.equal(length(filtered_fastq_files), 0))){
    stop("No filtered fastq files were found, problem in quality filtering process\n");
}
line_end <- length(filtered_fastq_files);


## Creating folder to perform the alignment

setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
system("mkdir alignment");
Sys.sleep(1);
setwd(paste0(fastq_files_folder, "/tmp/scripts/alignment"));


if(isTRUE(all.equal(species, "ant"))){ ## will use STAR aligner
    
    if(single_end){
        output_name_prefix <- unlist(strsplit(basename(filtered_fastq_files), ".fastq"))
        read_length_vec <- unlist(lapply(filtered_stat_files, function(x){
                            fi <- readLines(x)
                            fragLength <- unique(na.omit(as.numeric(strsplit(fi[grep("Maximum read length",fi)], " ")[[1]])))
                            return(fragLength)}))
        to_write_STARPaired <- paste(filtered_fastq_files, 1, 3, 3, 10, 1000000, output_name_prefix, 255, "None", read_length_vec, "Unstranded", sep=";")
        write(to_write_STARPaired, file="alignment.conf", ncolumns=1)
        command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/STAR_singleend_ant alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf");
    }else{
        
        fastq_factors <-as.factor(unlist(lapply(strsplit(filtered_fastq_files, "_R"), "[[", 1)))
        fastq_grouped_list <- split(filtered_fastq_files, fastq_factors)
        pasted_fastq_vec <- as.character(unlist(lapply(fastq_grouped_list, function(x){paste(x,collapse=";")})))
        output_name_prefix <- as.character(unlist(lapply(fastq_grouped_list, function(x){paste0(strsplit(basename(x), "_filtered")[[1]][1], "_R2")})))
        read_length_vec <- unlist(lapply(filtered_stat_files, function(x){
                    fi <- readLines(x)
                    fragLength <- unique(na.omit(as.numeric(strsplit(fi[grep("Maximum read length",fi)], " ")[[1]])))
                    return(fragLength)}))
        to_write_STARPaired <- paste(pasted_fastq_vec, 1, 3, 3, 10, 1000000, output_name_prefix, 255, "None", read_length_vec, "Unstranded", sep=";")
        write(to_write_STARPaired, file="alignment.conf", ncolumns=1)
        line_end <- line_end/2
        command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/STAR_pairedend_ant alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf");
    }
}else{
    if(chipseq){
        
        if(single_end){
            write(filtered_fastq_files, file="alignment.conf", ncolumns=1);
            
            if(isTRUE(all.equal(species, "mouse"))){
                command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_singleEnd_unireads_mouse alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf ", genome, " ", 3);
            }else{
                command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_singleEnd_unireads_human alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf ", genome, " ", 3);
            }
        }else{
            
            fastq_factors <-as.factor(unlist(lapply(strsplit(filtered_fastq_files, "_R"), "[[", 1)))
            fastq_grouped_list <- split(filtered_fastq_files, fastq_factors)
            to_write_bowtiePaired <- as.character(unlist(lapply(fastq_grouped_list, function(x){return(paste(x[1], x[2], sep=";"))})))
            write(to_write_bowtiePaired, file="alignment.conf", ncolumns=1);
            line_end <- line_end/2
            
            if(isTRUE(all.equal(species, "mouse"))){
                
                command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_pairedEnd_unireads_mouse alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf ", genome, " ", 3);
            }else{
                
                command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_pairedEnd_unireads_human alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf ", genome, " ", 3);
            }
        }
        
    }else{
        
        exp_name_vec <- unlist(lapply(strsplit(basename(filtered_fastq_files), "_filtered.fastq"),"[",1));
        write(paste(filtered_fastq_files, exp_name_vec, sep=";"), file="alignment.conf", ncolumns=1);
        
        if(isTRUE(all.equal(species, "mouse"))){  
            
            ##### !!!!!!!!!! TO MODIFY USING --read-edit-dist
            command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/tophat_SE_UR_submission_mouse alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " 2 ", " alignment.conf");
            
        }else{
            
            ##### !!!!!!!!!! TO MODIFY USING --read-edit-dist
            command_alignment <- paste0("/ifs/home/descon01/cluster/scripts/executable/tophat_SE_UR_submission_human alignment_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number," 2 ", " alignment.conf");
        }
    }
}


system(command_alignment);


if(chipseq && spikein){
    
    if(single_end)
        command_alignment_spikein <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_singleEnd_unireads_drosophila_forspikein alignment_", analysis_name, " alignment_spikein_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf dm3 3")
    else
        command_alignment_spikein <- paste0("/ifs/home/descon01/cluster/scripts/executable/bowtie_pairedEnd_unireads_drosophila_forspikein alignment_", analysis_name, " alignment_spikein_", analysis_name, " 1 ", line_end, " ", min_cpu_number, " ", max_cpu_number, " alignment.conf dm3 3")
    system(command_alignment_spikein);
}



##------------------ Launching the next step of pipeline

setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_2 <- "#!/bin/bash";
command_part_2 <- c(command_part_2, "#$ -cwd");
command_part_2 <- c(command_part_2, "module load r/3.2.0");
command_part_2 <- c(command_part_2, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part3.R --fastqFilesFolder ", fastq_files_folder, " --chipseq ", chipseq, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --fragmentSize ", fragment_size, " --wigFixed ", wig_fixed, " --wigVariable ", wig_variable, " --spikein ", spikein," --filteredFastqFiles ", paste(filtered_fastq_files, collapse=" "), " --analysisName ", analysis_name, " --minCpuNumber ", min_cpu_number, " --maxCpuNumber ", max_cpu_number));

write(command_part_2, file="load_part3.sh", ncolumns=1);

if(chipseq && spikein){
    command_part2 <- paste("qsub -q all.q -hold_jid alignment_spikein_", analysis_name, " load_part3.sh", sep="");
}else{
    command_part2 <- paste("qsub -q all.q -hold_jid alignment_", analysis_name, " load_part3.sh", sep="");
}

system(command_part2);

