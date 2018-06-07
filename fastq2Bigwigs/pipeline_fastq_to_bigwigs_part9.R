############
# This last part converts the spiked-in normalized exp to bigwigs.
# Descostes August 2017
############


library(Rargs);


################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


paramsDefinition[["--fastqFilesFolder"]] <- list(variableName="fastq_files_folder", numeric=F, mandatory=T, description="Folder paths to the fastq files.");
paramsDefinition[["--species"]] <- list(variableName="species", numeric=F, mandatory=T, description="Name of species. Currently human and mouse are supported.");

#fastq_files_folder <- "/ifs/home/descon01/data/data_august_2017/fasteq/orlando_spikein/0_percent_rep1";
#species <- "human"




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


wig_files <- list.files(paste0(fastq_files_folder, "/wig_files"), ".wig", full.names=TRUE);
chrom_size <- if(isTRUE(all.equal(species, "human"))) "/ifs/home/descon01/cluster/Annotations/human/hg19/hg19_chromsize_forbigwig.txt" else "/ifs/home/descon01/cluster/Annotations/mouse/mm10/chromsize_mm10_forbigwig.txt";

setwd(paste0(fastq_files_folder, "/tmp/scripts/"));
system("mkdir wig2bigwig");
Sys.sleep(3);
setwd(paste0(fastq_files_folder, "/tmp/scripts/wig2bigwig"));

write(paste(wig_files, chrom_size, sep=";"), file = "wig2bigwig.conf", ncolumns=1);
line_end <- length(wig_files);

command_wig2big <- paste0("/ifs/home/descon01/cluster/scripts/executable/wig2bigwig wig2bigwig 1 ", line_end, " wig2bigwig.conf");
system(command_wig2big);
