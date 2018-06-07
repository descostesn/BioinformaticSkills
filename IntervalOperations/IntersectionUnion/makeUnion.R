##############
# This script computes the union of intervals from several bed files and output a gff file.
# Descostes, September 2015
##############


library(IRanges);
library("Rargs");

################
# PARAMETERS
################


# Define the required options
paramsDefinition=list();

paramsDefinition[["--FilesPathVec"]]=list(variableName="files_path_vec", numeric=F, mandatory=T, description="A space separated list of bed or gff files to do the union on.");
paramsDefinition[["--outputPath"]]=list(variableName="output_path", numeric=F, mandatory=T, description="Single path to the output folder.");
paramsDefinition[["--outputFile"]]=list(variableName="output_file", numeric=F, mandatory=T, description="Single name of the output file.");
paramsDefinition[["--inputFormat"]]=list(variableName="input_format", numeric=F, mandatory=T, description="can be bed or gff.");


##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

if(input_format != "bed" && input_format != "gff")
{
    stop("input format should be bed or gff\n");
}


if(!file.exists(output_path))
{
    dir.create(output_path, recursive = TRUE)
}


# Reading bed files

cat("Reading files\n");

files_vec <- list();

for(i in 1:length(files_path_vec)) 
{
    files_vec[[i]] <- read.table(files_path_vec[i], stringsAsFactors=FALSE);
}

union_files <- do.call(rbind,files_vec);

if(input_format == "bed"){
    
    intervals_rangedData <- RangedData(IRanges(start = union_files[,2], end = union_files[,3]), space = union_files[,1]);
}else{
    intervals_rangedData <- RangedData(IRanges(start = union_files[,4], end = union_files[,5]), space = union_files[,1]);
}

cat("Reducing intervals\n");
intervals_rangedData <- reduce(intervals_rangedData);

gff_file_to_write <- cbind(seqName = as.character(intervals_rangedData[["space"]]), source = "union_files", feature = "union_files", start = start(intervals_rangedData), end = end(intervals_rangedData), 
        score = 0, strand = '+', frame = ".", group = make.unique(rep("group", length(start(intervals_rangedData))), sep="-"));

write.table(gff_file_to_write, file=paste(output_path, output_file, sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);



