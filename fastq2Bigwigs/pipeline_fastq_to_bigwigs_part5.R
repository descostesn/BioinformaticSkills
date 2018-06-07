#############
## This script loads the Pasha treatment of the sorted bam files.
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


## Preparing the conf file to load pasha

if(!isTRUE(all.equal(species, "ant"))){
    bam_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), "_SORTED_PICARD_COOR.bam", full.names = TRUE)
}else{
    bam_files_vec <- list.files(paste0(fastq_files_folder, "/bam_files"), "out.bam", full.names = TRUE)
}

if(isTRUE(all.equal(length(bam_files_vec), 0)))
    stop("Bam files not found, problem in the conversion of sam files.\n");

if(chipseq && spikein) ## Only performing the pipeline on endogenous data
{
    bam_files_vec <- bam_files_vec[-grep("dm3",bam_files_vec)];
}

if(chipseq)
{
    if(!isTRUE(all.equal(species, "ant")))
        expnames <- unlist(lapply(strsplit(basename(bam_files_vec), "_bwalign"), "[", 1))
    else
        expnames <- unlist(strsplit(basename(bam_files_vec), "Aligned.sortedByCoord.out.bam"))
}else
{
    expnames <- unlist(lapply(strsplit(basename(bam_files_vec), "_SORTED_"), "[", 1));
}

filepath <- bam_files_vec;
inputtype <- "BAM";

if(chipseq)
{
    arthefactth <- rep(7000000, length(filepath)); 
    elongationsize<- rep(fragment_size, length(filepath));
    
}else{
    arthefactth <- rep("NA", length(filepath));
    elongationsize<- rep(0, length(filepath));
}

pairedend<-!single_end;
subfolder<- "wig_files";
reportfolder<-"log_report";
wigfs <- wig_fixed;
wigvs<- wig_variable;
gff<-"FALSE";
nbcpu<- 1;
keeptmp<-"FALSE"
chr_prefix <- "chr";
annogff <- "NA";
multifile <- "NULL";
binsize <- 50;
midpoint <- "FALSE";

if(isTRUE(all.equal(species, "human"))){
    refile <- "/ifs/home/descon01/cluster/Annotations/human/hg19/hg19.ref";
}else if(isTRUE(all.equal(species, "mouse"))){
    refile <- "/ifs/home/descon01/cluster/Annotations/mouse/mm10/mm10.ref";
}else{
    refile <- "/ifs/home/descon01/cluster/Annotations/ants/Hsalv8_5/hsal85.ref"
}


to_write <- paste(expnames, filepath, inputtype, arthefactth, binsize, refile, pairedend, midpoint, subfolder, reportfolder, wigfs, wigvs, gff, elongationsize, nbcpu, keeptmp, chr_prefix, annogff, multifile, sep=";");

setwd(paste0(fastq_files_folder, "/tmp/scripts"));
system("mkdir pipeline_pasha");
Sys.sleep(5);
setwd(paste0(fastq_files_folder, "/tmp/scripts/pipeline_pasha"));

write(to_write, file = "pasha.conf", ncolumns = 1);
Sys.sleep(5);

line_end <- length(to_write);

command_pasha <- paste0("/ifs/home/descon01/cluster/scripts/executable/pipeline_pasha pipeline_pasha_", analysis_name, " 1 ", line_end, " pasha.conf");

system(command_pasha);



## ------------------------ Launch the next step


setwd(paste0(fastq_files_folder, "/tmp/scripts/steps_transition"));

command_part_5 <- "#!/bin/bash";
command_part_5 <- c(command_part_5, "#$ -cwd");
command_part_5 <- c(command_part_5, "module load r/3.2.0");
command_part_5 <- c(command_part_5, paste0("Rscript /ifs/home/descon01/cluster/scripts/R_scripts/core_facility_lab/fastq_to_bigwig/pipeline_fastq_to_bigwigs_part6.R --fastqFilesFolder ", fastq_files_folder, " --chipseq ", chipseq, " --genome ", genome, " --species ", species, " --singleEnd ", single_end, " --spikein ", spikein, " --analysisName ", analysis_name));

write(command_part_5, file="load_part6.sh", ncolumns=1);

command_part5 <- paste("qsub -q all.q -hold_jid pipeline_pasha_", analysis_name, " load_part6.sh", sep="");

system(command_part5);


