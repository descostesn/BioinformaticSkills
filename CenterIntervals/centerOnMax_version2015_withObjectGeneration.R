################
# This script combines the object generation (with save) and centering on max 
# values. It has 'version 2015' in the name because it does not use the new 
# rtracklayer package to extract values.
# Descostes March 2018
################


library(GenomicAnalysis)
library( "Rargs");

#################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--annotationFolder"]] <- list(variableName="annotationFolder", numeric=F, mandatory=T, description="Folder containing the annotations file.");
paramsDefinition[["--annotationFile"]] <- list(variableName="annotationFile", numeric=F, mandatory=T, description="Name of the annotations file (in the annotationFolder).");
paramsDefinition[["--annotationComment"]] <- list(variableName="annotationComment", numeric=F, mandatory=T, description="Comment on the type of annotation described in the annotationFile.");
paramsDefinition[["--annotationFormat"]] <- list(variableName="annotationFormat", numeric=F, mandatory=T, description="Format of the annotationFile ('bed', 'bedShort', 'gff'.");
paramsDefinition[["--wigFolders"]] <- list(variableName="wigFolder", numeric=F, mandatory=T, description="Folder containing the WIG file to treat.");
paramsDefinition[["--wigFiles"]] <- list(variableName="wigFile", numeric=F, mandatory=T, description="Unique WIG file name. The names are treated as regular expressions.");
paramsDefinition[["--expNames"]] <- list(variableName="expName", numeric=F, mandatory=T, description="Unique experiment name.");
paramsDefinition[["--annotationExtensions"]] <- list(variableName="annotationExtensions", numeric=T, mandatory=T, description="Unique extension applied at 5' and 3' of each (extended) annotation.");
paramsDefinition[["--spaceSizes"]] <- list(variableName="spaceSize", numeric=T, mandatory=T, description="Exclusion zone applied at 5' and 3' of each annotation.");
paramsDefinition[["--outputFolders"]] <- list(variableName="selectionFolder", numeric=F, mandatory=T, description="Unique folder where results will be output.");
paramsDefinition[["--binSize"]] <- list(variableName="binSize_value", numeric=T, mandatory=T, description="Value of the bin size of the wig files.");
paramsDefinition[["--nbCPU"]] <- list(variableName="nbCPU_value", numeric=T, mandatory=T, description="Number of CPU to use to generate the objects.");

# optional arguments
paramsDefinition[["--annotationFoldersExclude"]] <- list(variableName="annotationFoldersExclude", numeric=F, mandatory=F, description="Sapce separated list of folder containing annotation files to exclude", default= NA);
paramsDefinition[["--annotationFilesExclude"]] <- list(variableName="annotationFilesExclude", numeric=F, mandatory=F, description="Space separated list of annotation file names.  This list must be in the same order as the list provided as WigFoldersExclude.", default= NA);
paramsDefinition[["--annotationCommentExclude"]] <- list(variableName="annotationCommentExclude", numeric=F, mandatory=F, description="Space separated list of comment on the type of annotation to exclude This list must be in the same order as the list provided as WigFoldersExclude.", default= NA);
paramsDefinition[["--annotationFormatExclude"]] <- list(variableName="annotationFormatExclude", numeric=F, mandatory=F, description="Space separated lits of format of the annotation files This list must be in the same order as the list provided as WigFolders.", default= NA);
paramsDefinition[["--RNAseqFiles"]] <- list(variableName="RNAseqFileList", numeric=F, mandatory=F, description="List of files containing expression (rna-seq)", default=NA, postConversion=as.list);


#annotationFolder <- c("/ifs/home/descon01/analysis/K27M_project/overlap/peaks/data_feb2018/EZH2_8_timecourseK27M/peaks_per_circle/")
#annotationFile <- c("EZH2_0_EZH2_8_6H_EZH2_8_12H_EZH2_8_24H_EZH2_8_72H.gff")
#annotationComment <- c("commonEzh2inK27M")
#annotationFormat <- c("gff")
#wigFolder <- c("/ifs/home/descon01/data/data_february2018/bam_files/james_M_timecourse/Result_pasha_unireads_SE/AllReads/")
#wigFile <- c("WIGfs_EZH2_8_24hr-H3K27M-293TREX_filtered_unireads_elEst191_AThr2_bin50-RPM_BGSub-scaleReverse-Spikedin.bw")
#expName <- c("EZH2_8_24hr-H3K27M")
#annotationExtensions <- 15000
#spaceSize <- "NA"
#selectionFolder <- c("/ifs/home/descon01/analysis/K27M_project/objects/based_on_venndiagram/feb2018/EZH2_0_EZH2_8_6H_EZH2_8_12H_EZH2_8_24H_EZH2_8_72H_K27M/spNA/")
#binSize_value <- 50
#nbCPU_value <-10
#annotationFoldersExclude<-NA
#annotationFilesExclude <-NA
#annotationCommentExclude <- NA
#annotationFormatExclude<-NA
#RNAseqFileList<-NA
        



######################
# MAIN
######################

# Retreives the parameters
getParams(paramsDefinition);

if(!isTRUE(all.equal(length(annotationExtensions),1)) 
        || !isTRUE(all.equal(length(wigFolder), 1))
|| !isTRUE(all.equal(length(expName),1))
|| !isTRUE(all.equal(length(wigFile),1)))
    stop("This script treats only one experiment at a time")

if(!isTRUE(all.equal(length(spaceSize),1)))
    stop("Only one spacesize can be used")

if(is.na(RNAseqFileList))
    RNAseqFileList = list()


# Compose the annotation list
annotationsDescription <- list(folder=c(annotationFolder),
        file=c(annotationFile),
        annoType=c(annotationComment),
        fileType=c(annotationFormat))

# Compose the annotation to excludelist
annotationsToExclude <- list(folder= c(annotationFoldersExclude),
        file= c(annotationFilesExclude),
        annoType= c(annotationCommentExclude),	 
        fileType=c(annotationFormatExclude))


# Execute the Object generation on each entry
    
message("Generating the object")
        
selectionFile=paste("GenesSelect_", expName, "_spaceSize_", spaceSize, 
        "_ProLength_", annotationExtensions, sep="")

message("\t Selection File:", selectionFile)
message("\t SpaceSize", spaceSize)

current_object <- generate_genomic_information_object(expName=expName,
        wigFolder=wigFolder, 
        wigFile=wigFile, 
        annotationsDescription=annotationsDescription,
        annotationsToExclude=annotationsToExclude,
        exprsFiles=RNAseqFileList, 
        spaceSize=spaceSize, 
        profileLength=annotationExtensions,
        binSize=binSize_value,
        selectionOutputFolder=selectionFolder,
        selectionOutputFile=selectionFile,
        nbCPUs = nbCPU_value)
            
        

message("Centering on max value")

current_name <- expName
values_inside_list <- lapply(current_object, "[[", "probe.valueINSIDE")
start_list <- lapply(current_object, "[[", "gene.begin")
end_list <- lapply(current_object, "[[", "gene.end")
chr_list <- lapply(current_object, "[[", "chr")
strand_list <- lapply(current_object, "[[", "gene.strand")
name_list <- lapply(current_object, "[[", "gene.name")

multiple_max_count <- 0

max_index <- sapply(values_inside_list, function(x){
            
            index <- which(x == max(x));
            
            if(length(index) == 0)
            {
                stop("\n Problem in the data, no max was found\n\n");
            }
            else if(length(index) > 1)
            {
                cat("Several max were found, take a random one\n");
                multiple_max_count <<- multiple_max_count + 1;
                index <- index[sample(1:length(index),1)];
            }
            
            return(index);
        })

message("The number of annotations having more than one max is: ", 
        multiple_max_count, "/", length(chr_list))
start_vec <- unlist(start_list) + (max_index*binSize_value)
end_vec <- start_vec;

gff_df <- cbind(unlist(chr_list), paste(current_name, "-max_centered",sep=""),
        unlist(name_list), start_vec, end_vec, 0, unlist(strand_list), ".", 
        ".")

message("\t Writing output gff");
write.table(gff_df, file=paste(selectionFolder, "/", expName, 
                "-centeredmax.gff", sep=""), sep="\t", quote=F, row.names=F, 
        col.names=F);

