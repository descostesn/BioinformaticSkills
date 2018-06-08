########
# This script format the output of rsat dna pattern tool in gff format.
# IMPORTANT: lines beginning with a ";" should be removed before running this script.
# Descostes 27/02/2013
#######


### Includes
library(Rargs)


#################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


# Required arguments

paramsDefinition[["--folderVector"]] <- list(variableName="folderVector", numeric=F, mandatory=T, description="vector containing the folder path to the rsat files (space separated).");
paramsDefinition[["--annoFileVector"]] <- list(variableName="annoFileVector", numeric=F, mandatory=T, description="vector containing the file path of the semicolon separated annotations used for retrieving fasta sequences, these one were used with rsat to detect pattern (space separated).");
paramsDefinition[["--outputFolderVector"]] <- list(variableName="outputFolderVector", numeric=F, mandatory=T, description="vector containing the folders in which gff files will be written in (space separated).");
paramsDefinition[["--location"]] <- list(variableName="location", numeric=F, mandatory=T, description="Unique string telling where the pattern were detected. Should be TSS or TES.");
paramsDefinition[["--profileLengthBeforeVec"]] <- list(variableName="profileLengthBeforeVec", numeric=T, mandatory=T, description="vector of numeric values defining the size of the interval before the location used for detecting pattern with rsat. One numeric values per output folder should be given.");	
paramsDefinition[["--profileLengthAfterVec"]] <- list(variableName="profileLengthAfterVec", numeric=T, mandatory=T, description="vector of numeric values defining the size of the interval after the location used for detecting pattern with rsat. One numeric values per output folder should be given.");	



########
# MAIN
#######


# Retreives the parameters
getParams(paramsDefinition);


#Restrict the annotation of reference to the profile length

for(i in 1:length(folderVector))
{
	cat("Folder ", i, "/", length(folderVector),":\n");	
	motifFileVector <- list.files(folderVector[i], pattern=".txt", full.names = TRUE)
	annoFile <- annoFileVector[i];
	outputFolder <- outputFolderVector[i];
	profileLengthBefore <- profileLengthBeforeVec[i];
	profileLengthAfter <- profileLengthAfterVec[i];
	
	if(!file.exists(outputFolder))
	{
		dir.create(outputFolder, recursive=TRUE);
	}
	
	#Read the file of annotation and convert it to gff
	annoDataFrame <- read.table(annoFile, sep=":");
	annoDataFrame <- data.frame(seqname = annoDataFrame$V1, source = annoDataFrame$V4, feature = "conversion", start = annoDataFrame$V2, end = annoDataFrame$V3, score = 0, strand = annoDataFrame$V5, frame = '.', 
			group = make.unique(rep("group", length(annoDataFrame$V4)), sep="-"));
	
	#Restricting the annoDataFrame to the interval defined by the profileLength parameter
	
	for(i in 1:nrow(annoDataFrame))
	{
		if((as.character(annoDataFrame[i,"strand"]) == '-' && location == "TSS") || (as.character(annoDataFrame[i,"strand"]) == '+' && location == "TES"))
		{
			annoDataFrame[i,"start"] <- annoDataFrame[i,"end"];
		}
		else if((as.character(annoDataFrame[i,"strand"]) == '-' && location == "TES") || (as.character(annoDataFrame[i,"strand"]) == '+' && location == "TSS"))
		{
			annoDataFrame[i,"end"] <- annoDataFrame[i,"start"];
		}
		else
		{
			stop("\n Conditions are not respected for orientations, contact the developper\n")
		}
		
		if((annoDataFrame[i,"start"] - profileLengthBefore) < 0)
		{
			stop("\n Negative coordinates in the annotation file for ", annoDataFrame[i,],"\n\n");
		}
		
		annoDataFrame[i,"start"] <- annoDataFrame[i,"start"] - profileLengthBefore;
		annoDataFrame[i,"end"] <- annoDataFrame[i,"end"] + profileLengthAfter;
	}
	
	
	for(j in 1:length(motifFileVector))
	{
		cat("\t file", j, "/", length(motifFileVector), "\n");
		motifFilePath <- motifFileVector[j];
		elementPath <- strsplit(motifFilePath, split="/")[[1]];
		namefile <- strsplit(elementPath[length(elementPath)],split=".txt")[[1]];
		
		if(file.info(motifFilePath)$size != 0)
		{
			#read the motif file
			motifFile <- read.table(motifFilePath); #result of DNA pattern search
			
			#Retrieve the index of the corresponding annotation of annoDataFrame
			numberanno <- motifFile$V4;
			
			#Changing the motifFile dataFrame to gff format
			motifFile$V1 <- as.character(annoDataFrame[numberanno,1]) #replace the first column by the corresponding chrom
			motifFile$V2 <- as.character(annoDataFrame[numberanno,2]) #replace the second column by the annotation name
			motifFile$V3 <- "conversion";
			motifFile$V4 <- annoDataFrame[numberanno,4] + motifFile$V5 - 1; #put the begin of coordinates in position 4
			motifFile$V5 <- annoDataFrame[numberanno,4] + motifFile$V6 - 1; #put the end of coordinates in position 5
			motifFile$V6 <- 0;  #score of the gff file
			motifFile$V7 <- as.character(annoDataFrame[numberanno,7]); #put the strand of the annotation in position 7
			motifFile$V8 <- ".";
			motifFile <- cbind(motifFile, group = as.character(annoDataFrame[numberanno,9])); #add the group
			
			#Writting the gff table to the output folder
			cat("\t\t Writing table to output folder\n");
			write.table(motifFile, file = paste(outputFolder, namefile, ".gff", sep=""), quote=F, row.names=F, col.names=F);
		}
		else
		{
			cat("\t\t File is empty, no gff file created\n");
		}
	}
	
}	




##########################################################