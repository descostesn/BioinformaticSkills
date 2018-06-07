############
# This script performs differential alternative splicing by using DEXseq. The Bioconductor package DEXseq implements a method to test for differential exon usage in comparative RNA-Seq experiments.
# Descostes, Jan 2017
############

library("DEXSeq");
library("NGSprofiling");
library("BiocParallel");
library("Rargs");


################
# PARAMETERS
################



# Define the required options
paramsDefinition=list();

paramsDefinition[["--countFilesVec"]]=list(variableName="countFiles_vec", numeric=F, mandatory=T, description="A space separated vector of count files to compare.");
paramsDefinition[["--flattenedFile"]]=list(variableName="flattened_file", numeric=F, mandatory=T, description="GFF file of annotations obtained from GTF file in step 'perform_count_read'.");
paramsDefinition[["--rowNamesVec"]]=list(variableName="row_names_vec", numeric=F, mandatory=T, description="A vector giving the names of the samples.");
paramsDefinition[["--conditionVec"]]=list(variableName="condition_vec", numeric=F, mandatory=T, description="A vector giving the condition of each sample.");
paramsDefinition[["--libtypeVec"]]=list(variableName="libtype_vec", numeric=F, mandatory=T, description="A vector giving the library type which should be single-end or paired-end.");
paramsDefinition[["--outputFolder"]]=list(variableName="output_folder", numeric=F, mandatory=T, description="String giving the path to the output folder.");
paramsDefinition[["--fdrCutoff"]]=list(variableName="fdr_cutoff", numeric=T, mandatory=T, description="Single numeric giving the false discovery rate value. ex: 0.05");
paramsDefinition[["--nbCpu"]]=list(variableName="nb_cpu", numeric=T, mandatory=T, description="A single numeric giving the number of CPU.");


#countFiles_vec <- c("/usr/local/lib/R/site-library/pasilla/extdata/treated1fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/treated2fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/treated3fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/untreated1fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/untreated2fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/untreated3fb.txt", "/usr/local/lib/R/site-library/pasilla/extdata/untreated4fb.txt");
#flattened_file <- "/usr/local/lib/R/site-library/pasilla/extdata/Dmel.BDGP5.25.62.DEXSeq.chr.gff";
#row_names_vec <- c("treated1", "treated2", "treated3", "untreated1", "untreated2", "untreated3", "untreated4");
#condition_vec <- c("knockdown", "knockdown", "knockdown","control", "control", "control", "control");
#libtype_vec <- c("single-end", "paired-end", "paired-end", "single-end", "single-end", "paired-end", "paired-end"); 
#output_folder <- "/home/descostes/Documents/test/DEXseq/";
#fdr_cutoff <- 0.05;
#nb_cpu <- 2;


################

		
##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


checkingOutputFolder(output_folder);
BPPARAM <- MulticoreParam(workers=nb_cpu);

#Reading the data
countFiles <- countFiles_vec;

cat("The files selected for the analysis are:\n");
print(basename(countFiles));

flattenedFile <- flattened_file;
cat("The annotation file considered is:\n");
print(basename(flattenedFile));


#Building the sample table
sampleTable <- data.frame(row.names = row_names_vec, condition = condition_vec, libType = libtype_vec);


#Constructing the input object

cat("Loading the experiments\n");
dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= ~ sample + exon + condition:exon, flattenedfile=flattenedFile);


# Use the following commented lines to check the structure of your data
#split( seq_len(ncol(dxd)), colData(dxd)$exon )
#head( featureCounts(dxd), 5 )
#head( rowRanges(dxd), 3 )
#sampleAnnotation( dxd )


# Performing Normalization:

cat("Performing normalization of the data: Estimating size factors, dispersion, ");
#Different samples might be sequenced with different depths. In order to adjust for such coverage biases, we estimate size factors, which measure relative sequencing depth.
dxd <- estimateSizeFactors(dxd);

#To test for differential exon usage, we need to estimate the variability of the data. This is necessary to be able to distinguish technical and biological variation (noise) from real effects on exon usage due to the different conditions. The information on the strength of the noise is inferred from the biological replicates in the data set and characterized by the so-called dispersion.
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM);

png(filename=paste(output_folder, "dispersion.png",sep=""), width = 600, height = 600, bg = "transparent");
plotDispEsts(dxd);
dev.off();


# Estimating the differential exon usage
cat("Estimating the differential exon usage\n");
dxd <- testForDEU(dxd, BPPARAM=BPPARAM);
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM);


#Retrieving results
dxr1 <- DEXSeqResults(dxd);
nb_affected_exons <- table(dxr1$padj < fdr_cutoff)[2];
nb_affected_genes <- table(tapply(dxr1$padj < fdr_cutoff, dxr1$groupID, any))[2];
cat("The number of differentially used exons at FDR < ", fdr_cutoff*100, "%", " is ", nb_affected_exons, "\n");
cat("The number of genes affected at FDR < ", fdr_cutoff*100, "%", " is ", nb_affected_genes, "\n");

png(filename=paste(output_folder, "affected_exons.png",sep=""), width = 600, height = 600, bg = "transparent");
plotMA(dxr1, cex=0.8, main=paste(nb_affected_exons, " affected exons and ", nb_affected_genes, " affected genes \n", fdr_cutoff*100, " % FDR", sep=""));
dev.off();


# Plotting the differential exon usage for each gene
cat("Generating report in html format at: ", output_folder, "\n");
DEXSeqHTML(dxr1, color=c("#FF000080", "#0000FF80"), path=output_folder, file="differential_splicing.html", FDR= fdr_cutoff);
cat("\t done.\n");






