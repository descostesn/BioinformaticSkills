#########
# This script extract each category from results of rsat
#########

library(Rargs)


#################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


# Required arguments

paramsDefinition[["--workingDirVec"]] <- list(variableName="workingDir_vec", numeric=F, mandatory=T, description="Path to the working directory containing the patterns.");


#########
#MAIN
#########


# Retreives the parameters
getParams(paramsDefinition);

for(workingDir in workingDir_vec)
{
	cat("\n\n Analyzing ", workingDir, "\n\n");
	
	listFiles <- list.files(workingDir, recursive=TRUE, full.names=TRUE);
	
	listFiles <- listFiles[grep(listFiles, pattern="allpattern.save")];
	count <- 1;
	
	for(filepattern in listFiles)
	{
		cat("Reading pattern.save ", count, "/", length(listFiles), "\n");
		count <- count+1;
		
		filedata <- read.table(filepattern);
		motifVec <- as.character(unique(filedata$V1));
		outputfolder <- paste(dirname(filepattern),"/",sep="");
		
		countm <- 1;
		
		for(motif in motifVec)
		{
			cat("\t Reading motif", countm, "/", length(motifVec), "\n");
			countm <- countm + 1;
			
			filedatafiltered <- filedata[which(filedata$V1 == motif),];
			write.table(filedatafiltered, file=paste(outputfolder, motif, ".txt", sep=""), quote=F, row.names=F, col.names=F)
		}
	}
}

