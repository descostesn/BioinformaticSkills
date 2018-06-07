############
# This script converts wigs to bigwigs if no spikein. Otherwise load the spike-in scaling process.
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
paramsDefinition[["--singleEnd"]]=list(variableName="single_end", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiment is single or paired ended. Only single end is currently supported.", postConversion=as.logical);
paramsDefinition[["--spikein"]]=list(variableName="spikein", numeric=F, mandatory=T, description="Boolean (TRUE/FALSE) indicating if the experiments are spiked in.", postConversion=as.logical);
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="Name of the analysis. This will be used to create jobs names on the cluster.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
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


if(!spikein){
    
    wig_files <- list.files(paste0(fastq_files_folder, "/wig_files"), ".wig", full.names=TRUE);
    chrom_size <- if(isTRUE(all.equal(species, "human"))){"/ifs/home/descon01/cluster/Annotations/human/hg19/hg19_chromsize_forbigwig.txt"}else if(isTRUE(all.equal(species, "mouse"))){
                "/ifs/home/descon01/cluster/Annotations/mouse/mm10/chromsize_mm10_forbigwig.txt"}else "/ifs/home/descon01/cluster/Annotations/ants/Hsalv8_5/chomsize_hsal85_forbigwig.txt"
    
    setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
    system("mkdir wig2bigwig");
    Sys.sleep(3);
    setwd(paste0(fastq_files_folder, "/tmp/scripts/wig2bigwig"));
    
    write(paste(wig_files, chrom_size, sep=";"), file = "wig2bigwig.conf", ncolumns=1);
    line_end <- length(wig_files);
    
    command_wig2big <- paste0("/ifs/home/descon01/cluster/scripts/executable/wig2bigwig wig2bigwig_", analysis_name, " 1 ", line_end, " wig2bigwig.conf");
    system(command_wig2big);
    
}else{
    
    bam_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), "_SORTED_PICARD_COOR.bam$", full.names = TRUE);
    indexes_dm3 <- grep("dm3", bam_files_vec);
    bam_files_endo <- bam_files_vec[-indexes_dm3];
    bam_files_exo <- bam_files_vec[indexes_dm3];
    
    bai_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), "bai$", full.names = TRUE);
    indexes_dm3 <- grep("dm3", bai_files_vec);
    bai_files_endo <- bai_files_vec[-indexes_dm3];
    bai_files_exo <- bai_files_vec[indexes_dm3];
    
    expnames_endo <- paste(unlist(lapply(strsplit(basename(bai_files_endo),"_bwalign"),"[",1)), genome,sep="_");
    expnames_exo <- paste(unlist(lapply(strsplit(basename(bai_files_exo),"_bwalign"),"[",1)), "dm3",sep="_");
    
    outputfile_endo <- paste("scores_", genome, ".txt", sep="");
    outputfile_exo <- "scores_dm3.txt";

    setwd(paste0(fastq_files_folder, "/tmp/scripts/spikeIn_process"));
    
    param_conf_file_endo <- paste(paste(bam_files_endo, collapse=" "), paste(bai_files_endo, collapse=" "), paste(expnames_endo, collapse=" "), outputfile_endo, if(single_end) "FALSE" else "TRUE", sep=";");
    param_conf_file_exo <- paste(paste(bam_files_exo, collapse=" "), paste(bai_files_exo, collapse=" "), paste(expnames_exo, collapse=" "), outputfile_exo, if(single_end) "FALSE" else "TRUE", sep=";");
    write(param_conf_file_endo, file="endo.conf", ncolumns=1);
    write(param_conf_file_exo, file="exo.conf", ncolumns=1);
    
    line_end <- 1; ## Since endo.conf and exo.conf will always have one line in the context of this pipeline
    
    command_endo <- paste0("/ifs/home/descon01/cluster/scripts/executable/compute_scaling_factors SCFendo_", analysis_name, " 1 ", line_end, " endo.conf");
    command_exo <- paste0("/ifs/home/descon01/cluster/scripts/executable/compute_scaling_factors_withhold SCFendo_", analysis_name, " SCFexo_", analysis_name, " 1 ", line_end, " exo.conf");
    
    system(command_endo);
    system(command_exo);
    
    
    #------------- launch next script
    
    setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));
    
    command_part_7 <- "#!/bin/bash";
    command_part_7 <- c(command_part_7, "#$ -cwd");
    command_part_7 <- c(command_part_7, "module load r/3.2.0");
    command_part_7 <- c(command_part_7, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part8.R --fastqFilesFolder ", fastq_files_folder, " --genome ", genome, " --species ", species, " --analysisName ", analysis_name));
    
    write(command_part_7, file="load_part8.sh", ncolumns=1);
    
    command_part7 <- paste("qsub -q all.q -hold_jid SCFexo_", analysis_name, " load_part8.sh", sep="");
    system(command_part7);
}



