library(Pasha)

library(Rargs)	


### Parameters from command line using RIO ###
# Parameters format definition
paramsDefinition=list()

# mandatory (MUST BE SPECIFIED BY USER, DEFAULT VALUES ARE SPECIFIED FOR INFORMATION)
paramsDefinition[["--expName"]]=list(variableName="expName", numeric=F, mandatory=T, description="Experiment name (used for output filenames)");
paramsDefinition[["--inputFile"]]=list(variableName="inputFile", numeric=F, mandatory=T, description="Aligned reads file name");
paramsDefinition[["--inputType"]]=list(variableName="fileType", numeric=F, mandatory=T, description="Aligned reads file format", default="BAM");
paramsDefinition[["--threshold"]]=list(variableName="artefactThr", numeric=T, mandatory=T, description="Threshold for artefact removal (NA for no removal, multiple values allowed, negative value for direct specification of repeats number allowed)", default=7000000);
paramsDefinition[["--bins"]]=list(variableName="binSize", numeric=T, mandatory=T, description="Size of the bins (for wig and GFF output, multiple values allowed)", default=50);
paramsDefinition[["--annotationGenomeFiles"]]=list(variableName="annotation_genome_files", numeric=F, mandatory=T, description=" A single path to a genome reference file or the ID of a genome for which a reference is bundled in the package (see details)");


# non mandatory with default, "expert mode"
paramsDefinition[["--pairedEnds"]]=list(variableName="pairedEndsExperiment", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), should the experiment be considered as a paired-ends experiment", default="FALSE", postConversion=as.logical);
paramsDefinition[["--midPoint"]]=list(variableName="midPoint", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), define the type of piling (midpoint or normal)", default="FALSE", postConversion=as.logical);
paramsDefinition[["--resultSubFolder"]]=list(variableName="resultSubFolder", numeric=F, mandatory=F, description="Subfolders containing results in the experiment folder", default="ResultPasha");
paramsDefinition[["--reportFilesSubFolder"]]=list(variableName="reportFilesSubFolder", numeric=F, mandatory=F, description="Subfolder containing Logs in the Result folder", default="ReportFiles");
paramsDefinition[["--WIGfs"]]=list(variableName="WIGfs", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), write a WIG fixed step output file", default="TRUE", postConversion=as.logical);
paramsDefinition[["--WIGvs"]]=list(variableName="WIGvs", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), write a WIG variable step output file", default="TRUE", postConversion=as.logical);
paramsDefinition[["--GFF"]]=list(variableName="GFF", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), write a GFF output file", default="FALSE", postConversion=as.logical);
paramsDefinition[["--elongationSize"]]=list(variableName="elongationSize", numeric=T, mandatory=F, description="Fragments size for tag elongation (NA for automatic estimation, 0 for no elongation)", default=NA);
paramsDefinition[["--elongationRangeMin"]]=list(variableName="elongationRangeMin", numeric=T, mandatory=F, description="Minimal value for estimation of fragments size", default=10);
paramsDefinition[["--elongationRangeMax"]]=list(variableName="elongationRangeMax", numeric=T, mandatory=F, description="Maximal value for estimation of fragments size", default=450);
paramsDefinition[["--elongationResolution"]]=list(variableName="elongationResolution", numeric=T, mandatory=F, description="Resolution for estimation of fragments size", default=10);
paramsDefinition[["--ignoreChr"]]=list(variableName="removeChrNamesContaining", numeric=F, mandatory=F, description="Regular expression identifying undesired chromosomes in output", default="random|hap");
paramsDefinition[["--ignoreInsertsOver"]]=list(variableName="ignoreInsertsOver", numeric=T, mandatory=F, description="In case of paired-ends experiments, specify the max insert size to be considered as valid", default=500);
paramsDefinition[["--CPU"]]=list(variableName="nbCPU", numeric=T, mandatory=F, description="Number of threads to start (AT THE COST OF MEMORY)", default=1);
paramsDefinition[["--keepTemp"]]=list(variableName="keepTemp", numeric=F, mandatory=F, description="Boolean (TRUE/FALSE), keep the intermediate results (rehabilitation steps)", default=TRUE, postConversion=as.logical);
paramsDefinition[["--inputChrPrefix"]]=list(variableName="inputChrPrefix", numeric=F, mandatory=F, description="Prefix used for chromosome number in reference genome", default="chr");
paramsDefinition[["--inputChrSuffix"]]=list(variableName="inputChrSuffix", numeric=F, mandatory=F, description="Suffix used for chromosome number in reference genome", default="''");
paramsDefinition[["--annotationFilesGFF"]] <- list(variableName="annotation_file_gff", numeric= F, mandatory=F, description= "A complex parameter (see details). A named vector of gff file paths. If this argument and annotationGenomeFiles are provided, the pipeline will generate or each rangeSelection (and total) a pdf file summarizing reads occupancy among annotations in gff files.", default= NA);

# multiread options
multireadSignalFile = NULL
paramsDefinition[["--multireadSignalFile"]]=list(variableName="multireadSignalFile", numeric=F, mandatory=F, description="Path to the file containing the distribution of the signal of multiread tags.");

# Retreives the parameters
getParams(paramsDefinition);

### Check that input file exists
if(!file.exists(inputFile)) stop("Input file does not exist !");

### Build a list for the experiment
experiment=list();
experiment[[expName]]=list(folderName=paste(dirname(inputFile),"/", sep=""), fileName=basename(inputFile), fileType=fileType, chrPrefix=inputChrPrefix, chrSuffix=inputChrSuffix, pairedEnds=pairedEndsExperiment, midPoint=midPoint)

### Build a list for the multiread signal file
multireadSignalFileList=list()
if( !is.null( multireadSignalFile))
{
	multireadSignalFileList[[expName]]=multireadSignalFile
}

# Do it for 50bp in priority
returnBinary=processPipeline(
		# I/O GENERAL PARAMETERS
		INPUTFilesList					= experiment,
		resultSubFolder					= resultSubFolder,
		reportFilesSubFolder			= reportFilesSubFolder,
		WIGfs							= WIGfs,
		WIGvs							= WIGvs,
		GFF								= GFF,
		# COMPLEX PARAMETERS (ATOMIC OR VECTORS OR LIST OF IT)
		incrArtefactThrEvery			= artefactThr,
		binSize							= binSize,
		elongationSize					= elongationSize,
		rangeSelection					= IRanges(0,-1),
		# ATOMIC PARAMETERS
		elongationEstimationRange		= c(mini=elongationRangeMin, maxi=elongationRangeMax, by=elongationResolution),
		rehabilitationStep				= c("orphans","orphansFromArtefacts"),
		removeChrNamesContaining		= removeChrNamesContaining,
		ignoreInsertsOver				= ignoreInsertsOver,
		nbCPUs							= nbCPU,
		keepTemp						= keepTemp,
		logTofile						= NULL,
		eraseLog						= TRUE,
		multiLocFilesList				= multireadSignalFileList,
		annotationFilesGFF				= annotation_file_gff,
		annotationGenomeFiles			= annotation_genome_files);



