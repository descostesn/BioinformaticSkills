###############
# This script converts a wig file in variable step to fixed step.
# Descostes September 2015
###############

library( "Rargs");

#################
# PARAMETERS
################


paramsDefinition <- list();

paramsDefinition[["--wigFolder"]] <- list(variableName="wigFolder", numeric=F, mandatory=T, description="Folder containing the wig file to convert");
paramsDefinition[["--wigFile"]] <- list(variableName="wigFile", numeric=F, mandatory=T, description="Wig file to convert");
paramsDefinition[["--nameWig"]] <- list(variableName="nameWig", numeric=F, mandatory=T, description="Name to write in the future converted wig file");
paramsDefinition[["--binning"]] <- list(variableName="binning", numeric=T, mandatory=T, description="Size of bins of the future converted wig file");
paramsDefinition[["--refFile"]] <- list(variableName="refile", numeric=F, mandatory=T, description="File path containg genome information (nb of chrom, names of chrom, length of chrom)");


#wigFolder <- "/home/descostes/Documents/analysis/EED-mutant/wig_files/";
#wigFile <- "WIGvs_E14_Y358A_H3K27me3_dm3_unireads_elManual156_AThr1_countsNorm.wig";
#nameWig <- "E14_Y358A_H3K27me3_dm3_unireads_elManual156_AThr1_countsNorm";
#binning <- 50;
#refile <- "/home/descostes/Documents/Annotations/drosophila/dm3/dm3.ref";







################
# FUNCTION
################


extract_start_end_vector <- function(last_position_value, length_current_chrom, position_vec, binding_values_vec)
{
	if(last_position_value <= length_current_chrom)
	{
		if(last_position_value == length_current_chrom)
		{
			position_star_vec <- position_vec[-length(position_vec)];
			position_end_vec <- position_vec[-1];
			corrected_bindingValues <- binding_values_vec[-length(position_vec)];
		}
		else
		{
			position_star_vec <- position_vec;
			position_end_vec <- c(position_vec[-1], length_current_chrom);
			corrected_bindingValues <- binding_values_vec;
		}
		
	}
	else
	{
		stop("Some values are referenced after the end of the chromosome, problem in the construction of your wig file\n\n");
	}
	
	return(list(position_star_vec, position_end_vec, corrected_bindingValues));
}


################



##########
# MAIN
##########

# Retreives the parameters
getParams(paramsDefinition);


#Reading the wig file
cat("Reading wig... \n")
wigFile <- readLines(paste(wigFolder, wigFile, sep=""));
count <- 1;
fixedList <- list();


#Reading the reference file

ref_file <- readLines(refile);
chrom_names <- unlist(strsplit(ref_file[3], " "));
chrom_length <- as.numeric(unlist(strsplit(ref_file[2], " ")));

cat("Retrieving the track line and variableStep lines...\n");

indexTrack <- grep("track", wigFile);
indexVariableStepLines <-  grep("variableStep", wigFile);

insert.at <- function(a, pos, values){
	dots <- values
	stopifnot(length(dots)==length(pos))
	result <- vector("list",2*length(pos))
	result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% pos))
	result[c(FALSE,TRUE)] <- dots
	unlist(result)
}



if(length(indexTrack) == 0) #If the line "track" is absent before the line "VariableStep", it is inserted. Both indexes are computed again after this operation 
{
	cat("Adding a track line\n");
	pos <- indexVariableStepLines;
	values <- rep("track type=wiggle_0", length(pos));
	wigFile <- insert.at(wigFile, pos,values);
	wigFile <- c("track type=wiggle_0", wigFile);
	wigFile <- wigFile[-length(wigFile)];
	
	indexTrack <- grep("track", wigFile);
	indexVariableStepLines <-  grep("variableStep", wigFile);
}

#changing the track character to adapt it to fixed format
trackLine <- paste("track type=wiggle_0 name=\"",  nameWig, "\"", sep="");

#changing the variableStep line to fit the fixedStep format
fixedStepLines <- sub("variableStep", "fixedStep", wigFile[indexVariableStepLines]);
fixedStepLines <- paste(fixedStepLines, " start=1 step=", binning, sep="");

#Retrieving the names of chrom in wig file
chrom_names_wig <- unlist(lapply(strsplit(unlist(lapply(strsplit(fixedStepLines,"chrom="),"[[",2))," start"),"[[",1));

cat("Retrieving the values and positions and transforming into fixed step...\n");

##Create interval of values
indexStartInterval <- indexVariableStepLines[1:(length(indexVariableStepLines)-1)] +1;
indexEndInterval <- indexTrack[2:length(indexTrack)] - 1;

# Adding the last interval
indexStartInterval <- c(indexStartInterval, indexVariableStepLines[length(indexVariableStepLines)]+1);
indexEndInterval <- c(indexEndInterval, length(wigFile));

matrixInterval <- cbind(indexStartInterval, indexEndInterval);

##Retrieve values in these intervals ie the values between the variableStep lines

#valuesList <- list()
valuesList <- apply(matrixInterval, MARGIN=1, function(x){ return(wigFile[x[1]:x[2]]);});

values_binned <- list();

##changing the values in variable steps into fixed steps

for(i in 1:length(valuesList)) 
{
	#split each element of x into a list. each element of the new list contain an index and a value
	current_chrom_name <- chrom_names_wig[i];
	cat("Treating chromosome: ", current_chrom_name, "\n");
	splittedIndexAndValue <- strsplit(valuesList[[i]], "\t");
	position_vec <- as.numeric(unlist(lapply(splittedIndexAndValue, "[[", 1)));
	binding_values_vec <- as.numeric(unlist(lapply(splittedIndexAndValue, "[[", 2)));
	index_ref_chrom_name <- which(chrom_names == current_chrom_name);
	length_current_chrom <- chrom_length[index_ref_chrom_name];
	last_position_value <- position_vec[length(position_vec)];
	
	if(position_vec[1] > 1)
	{
		position_vec <- c(1, position_vec);
		binding_values_vec <- c(0, binding_values_vec);
	}
	else if(position_vec[1] == 0)
	{
		stop("\n The first position is equal to zero, pb in the wig file\n");	
	}
	
	
	result <- extract_start_end_vector(last_position_value, length_current_chrom, position_vec, binding_values_vec);
	position_star_vec <- result[[1]];
	position_end_vec <- result[[2]];
	binding_values_vec <- result[[3]];
	
	length_vec <- position_end_vec - position_star_vec;
	
	if(length(length_vec) != length(binding_values_vec))
	{
		stop("\n The number of intervals should be equal to the number of values, pb in script, contact the developper\n");
	}
	
	valuesRepeated_vector <- unlist(mapply(function(x,y){return(rep(x,y));}, binding_values_vec, length_vec));
	
	#performing the bining of the values
	steps <- seq(1,length(valuesRepeated_vector),binning);
	
	values_binned[[i]] <- sapply(steps, function(step,orig)
			{
				nextStep <- step+binning-1;
				if(nextStep > length(orig))
				{
					nextStep <- length(orig)
				}
				
				return(mean(orig[step:nextStep]))
			},valuesRepeated_vector)
	
}

			
## Computing the list to write in the output file
fixedVector <- vector();

if(length(values_binned) != length(fixedStepLines))
{
	stop("\n The number of blocks of values should be equal to the number of lines containing the key word 'variableStep', problem in your wig file\n");
}

cat("Organizing all information before writting the file...\n");

for(i in 1:length(values_binned))
{
	fixedVector <- append(fixedVector, trackLine);
	fixedVector <- append(fixedVector, fixedStepLines[i]);
	fixedVector <- append(fixedVector, unlist(values_binned[[i]]));
}


outputFile <- paste(wigFolder, "WIGfs_", nameWig, "_bin", binning, ".wig", sep="")

cat("Writting wig... \n");
write(fixedVector, file= outputFile, ncolumns=1);
