#############################
# This script creates a matrix of colored nucleotides according to coordinates found with Rsat tools. It needs the file of motif coordinates and the file of semicolon sorted genes.
# Descostes, ap 2016
#############################

library(NGSprofiling);
library(GenomicRanges);
library(Rargs);

paramsDefinition <- list();


# Required arguments
paramsDefinition[["--annoFile"]] <- list(variableName="annoFile_path", numeric=F, mandatory=T, description="File of semi column separated annotations");
paramsDefinition[["--coordinatesForColoring"]] <- list(variableName="coordinatesForColoring", numeric=F, mandatory=T, description="Gff file containing the intervals to color", default=NA);
paramsDefinition[["--profileLengthBefore"]] <- list(variableName="profileLengthBefore", numeric=T, mandatory=T, description="coordinate of profile before the location defined");
paramsDefinition[["--profileLengthAfter"]] <- list(variableName="profileLengthAfter", numeric=T, mandatory=T, description="coordinate of profile after the location defined");
paramsDefinition[["--location"]] <- list(variableName="location", numeric=F, mandatory=T, description="Should be TSS");
paramsDefinition[["--outputFolder"]] <- list(variableName="outputFolder", numeric=F, mandatory=T, description="Folder in which the treeview file will be written");
paramsDefinition[["--nameExp"]] <- list(variableName="nameExp", numeric=F, mandatory=T, description="name of the experiment, generally the name of the gff file used");


################





##############
# MAIN
##############

getParams(paramsDefinition);

checkingOutputFolder(outputFolder);

if(length(annoFile_path) != 1)
{
	stop("\n The annotation file of the sorted genes should be unique\n");
}

if(location != "TSS")
{
	stop("\n This script only deals with TSS\n");
}

cat("Reading the file of sorted genes\n");
annoFile <- read.table(annoFile_path, stringsAsFactors=FALSE, sep=":");

cat("Transforming coordinates to TSS-profileLengthBefore and TSS+profileLengthAfter\n");
annoFile$V3 <- annoFile$V2 + profileLengthAfter;
annoFile$V2 <- annoFile$V2 - profileLengthBefore;

length_seq <- unique(annoFile$V3 - annoFile$V2); 

if(length(length_seq) != 1)
{
	stop("\n All sequences submitted should have the same length to draw the heatmap\n");
}


annoFile_rangedData <- gffToRangedData(data.frame(seqname=annoFile$V1, source="color_gff", feature=annoFile$V4,start=annoFile$V2,end=annoFile$V3,score=0,strand=annoFile$V5,frame=".", group="."))

for(i in 1:length(coordinatesForColoring)) 
{
	cat("Reading the file of coordinates", i, "/", length(coordinatesForColoring), "\n");
	current_coord <- read.table(coordinatesForColoring[i], stringsAsFactors=FALSE);
	current_coord_rangedData <- gffToRangedData(data.frame(seqname=current_coord$V1, source="color_gff", feature=current_coord$V2,start=current_coord$V4,end=current_coord$V5,score=0,strand=current_coord$V7,frame=".", group="."));
	
	cat("Retrieving the coordinates corresponding to each annotation by overlap\n");
	hits_result <- findOverlaps(current_coord_rangedData, annoFile_rangedData);
	
	#Separating each motif coordinate according to overlap with reference sequences
	list_coordinates <- split(queryHits(hits_result), subjectHits(hits_result));
	
	#Retrieving the start and end of all coordinates
	start_coord_vec <- start(current_coord_rangedData);
	end_coord_vec <- end(current_coord_rangedData);
	start_ref_seq <- start(annoFile_rangedData);
	
	count_ref_seq <- 1;  #the split follows the number of the factor, 1,2,...., so an increasing index is used to retrieve the start of the ref sequence
	
	cat("Building the matrix to generate heatmaps\n");
	matrix_list_result <- lapply(list_coordinates, function(x, start_vec, end_vec, start_ref, seq_length){
				
				start_reference <- start_ref[count_ref_seq];
				current_start <- start_vec[x]-start_reference;
				current_end <- end_vec[x]-start_reference;
				
				#Generating all indexes to replace by 1
				index_1 <- unlist(mapply(function(x,y){return(seq(x,y))}, current_start, current_end));
				
				result <- rep(0, seq_length);
				
				if(max(index_1) > seq_length) #This should not happen
				{
					stop("\n pb in the script\n");
				}
				
				result[index_1] <- 1;
				count_ref_seq <<- count_ref_seq + 1;
				return(result);
				
			}, start_coord_vec, end_coord_vec, start_ref_seq, length_seq);
	
	result_matrix <- do.call(rbind, matrix_list_result);
	
	cat("Writting result matrix \n");
	write.table(result_matrix, file= paste(outputFolder, "matrix-", nameExp[i], ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE);
	
	cat("Writting the files for visualization in treeview...");
	interpolatedCoordinates <- 1:ncol(result_matrix);
	annoNames <- read.table(annoFile_path); #Read the semicolon separated file and use each line as row name for the treeview file
	writeTreeviewCluster(paste(outputFolder, nameExp[i], sep=""), result_matrix, interpolatedCoordinates, annoNames, compress=FALSE);
}






