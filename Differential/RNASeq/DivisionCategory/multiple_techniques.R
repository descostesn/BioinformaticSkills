###########
# This script plots the distribution of gene expression from a bam file. define different classes with different classification techniques, and output bed files for each category.
# Descostes feb 2017
###########


library("NGSprofiling");
library("GenomicAlignments");
library("GenomicFeatures");
library("edgeR");
library("DESeq");
library("classInt");
library("Rargs");



################
# PARAMETERS
################


paramsDefinition <- list();


# Required arguments
paramsDefinition[["--bamFile"]] <- list(variableName="bam_file", numeric=F, mandatory=T, description="Single bam file from which the different categories of expression will be extracted.");
paramsDefinition[["--gtfFile"]] <- list(variableName="gtf_file", numeric=F, mandatory=T, description="GTF file used to compute reads count.");
paramsDefinition[["--singleEnd"]] <- list(variableName="single_end", numeric=F, mandatory=T, description="Boolean indicating if bam file contains sine or paired end tags", postConversion=as.logical);
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Vector of output folders corresponding to each set of bam files.");
paramsDefinition[["--bedAnno"]] <- list(variableName="bed_anno", numeric=F, mandatory=T, description="Bed file of annotations corresponding to the gtf file.");
paramsDefinition[["--nbOfGroups"]] <- list(variableName="nb_of_groups", numeric=T, mandatory=T, description="Nb of groups for classifying data.");



#bam_file <- c("/ifs/home/descon01/data/data_february2017/bam_files/gary_rnaseq_293/293TRex_rep1tech1_SORTED_PICARD_COOR.bam");
#gtf_file <- "/ifs/home/descon01/cluster/Annotations/human/hg19/gtf_files/refseq_hg19.gtf";
#single_end <- TRUE;
#output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/expression_distribution/rep1tech1/";
#bed_anno <- "/ifs/home/descon01/cluster/Annotations/human/hg19/bed_files/refseq_hg19.bed"; 
#nb_of_groups <- 3;

################




################
# FUNCTION
################


plot_histogram_table <- function(output_folder, file_name, data_mat, norm_name)
{
	png(filename=paste(output_folder, file_name, ".png",sep=""), width = 1000, height = 1000, bg = "transparent")
	hist(data_mat, xlab="log2(Expression)", ylab="nb of genes", breaks=max(data_mat), main=norm_name);
	dev.off();
}


highlight_breaks_and_nb_indiv <- function(class_interval, data_mat, output_folder, file_name)
{
	#Highlight the breaks on the histogram and plot the number of genes per categories
	
	palette_vec <- attr(findColours(class_interval,pal1),"palette");
	colors_vec <- vector();
	
	for(i in 1:(length(class_interval$brks)-1))
	{
		colors_vec <- c(colors_vec, rep(palette_vec[i], (class_interval$brks[i+1]-class_interval$brks[i])));
	}
	
	tmp_breaks <- hist(data_mat, breaks=max(data_mat),plot=FALSE);
	tmp_breaks <- sum(abs(min(tmp_breaks$mids)), abs(max(tmp_breaks$mids)));
	
	png(filename=paste(output_folder, file_name,".png",sep=""), width = 600, height = 600, bg = "transparent");
	hist(data_mat, xlab="log2(Expression)", ylab="nb of genes", breaks=tmp_breaks, main="RPKM", col=colors_vec);
	abline(v=class_interval$brks, col="blue");
	dev.off();
	
	png(filename=paste(output_folder, file_name,"-nb_elements.png",sep=""), width = 600, height = 600, bg = "transparent");
	barplot(print(class_interval), main= file_name, ylab="Nb of genes", col=palette_vec);
	dev.off();
}



separate_bed_by_class <- function(classInterval, data_mat, bed_fi, output_folder, name_classInterval)
{
	for(i in 1:(length(classInterval$brks)-1)) 
	{
		index_current_class <- which(data_mat > classInterval$brks[i] & data_mat <= classInterval$brks[i+1]);
		bed_current_class <- bed_fi[index_current_class,];
		write.table(bed_current_class, file=paste(output_folder, name_classInterval, "-category", i, ".bed", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE);
	}
}



computing_class_intervals <- function(data_mat, nb_of_groups, output_folder, bed_fi)
{
	cat("\t Computing data classes\n");
	# Computing different class definition
	#jenks_result <- classIntervals(edgeR_rpkm_matrix, n=nb_of_groups, style="jenks", dataPrecision=0, unique=TRUE);
	sd_result <- classIntervals(data_mat, n=nb_of_groups, style="sd", dataPrecision=0, unique=TRUE);
	equal_result <- classIntervals(data_mat, n=nb_of_groups, style="equal", dataPrecision=0, unique=TRUE);
	pretty_result <- classIntervals(data_mat, n=nb_of_groups, style="pretty", dataPrecision=0, unique=TRUE);
	quantile_result <- classIntervals(data_mat, n=nb_of_groups, style="quantile", dataPrecision=0, unique=TRUE);
	kmean_result <- classIntervals(data_mat, n=nb_of_groups, style="kmeans", dataPrecision=0);
	hclust_result <- classIntervals(data_mat, n=nb_of_groups, style="hclust", method="single", dataPrecision=0);
	fisher_result <- classIntervals(data_mat, n=nb_of_groups, style="fisher", dataPrecision=0, unique=TRUE);
	
	cat("\t Creating output folders\n");
	
	output_folder_sd <- paste(output_folder, "sd/", sep="");
	output_folder_equal <- paste(output_folder, "equal/", sep="");
	output_folder_pretty <- paste(output_folder, "pretty/", sep="");
	output_folder_quantile <- paste(output_folder, "quantile/", sep="");
	output_folder_kmean <- paste(output_folder, "kmean/", sep="");
	output_folder_hclust <- paste(output_folder, "hclust/", sep="");
	output_folder_fisher <- paste(output_folder, "fisher/", sep="");
	checkingOutputFolder(output_folder_sd);
	checkingOutputFolder(output_folder_equal);
	checkingOutputFolder(output_folder_pretty);
	checkingOutputFolder(output_folder_quantile);
	checkingOutputFolder(output_folder_kmean);
	checkingOutputFolder(output_folder_hclust);
	checkingOutputFolder(output_folder_fisher);
	
	cat("\t Plotting classes read out\n");
	
	highlight_breaks_and_nb_indiv(sd_result, data_mat, output_folder_sd, "sd");
	highlight_breaks_and_nb_indiv(equal_result, data_mat, output_folder_equal, "equal");
	highlight_breaks_and_nb_indiv(pretty_result, data_mat, output_folder_pretty, "pretty");
	highlight_breaks_and_nb_indiv(quantile_result, data_mat, output_folder_quantile, "quantile");
	highlight_breaks_and_nb_indiv(kmean_result, data_mat, output_folder_kmean, "kmean");
	highlight_breaks_and_nb_indiv(hclust_result, data_mat, output_folder_hclust, "hclust");
	highlight_breaks_and_nb_indiv(fisher_result, data_mat, output_folder_fisher, "fisher");
	
	
	cat("\t Synchronizing the bed file of annotation on the matrix\n");
	
	names_count_table <- rownames(data_mat);
	indexes_for_ordering <- match(names_count_table, bed_fi$V4);
	if(length(which(is.na(indexes_for_ordering))) != 0)
	{
		stop("\n The gtf and bed files provided contain different annotations. Make sure they do.\n");
	}
	if(length(indexes_for_ordering) != length(names_count_table))
	{
		stop("\n The names of annotations should be unique in the bed and gtf files\n");
	}
	bed_fi <- bed_fi[indexes_for_ordering,];
	
	
	cat("\t Separating genes by expression classes\n");
	
	separate_bed_by_class(sd_result, data_mat, bed_fi, output_folder_sd, "sd");
	separate_bed_by_class(equal_result, data_mat, bed_fi, output_folder_equal, "equal");
	separate_bed_by_class(pretty_result, data_mat, bed_fi, output_folder_pretty, "pretty");
	separate_bed_by_class(quantile_result, data_mat, bed_fi, output_folder_quantile, "quantile");
	separate_bed_by_class(kmean_result, data_mat, bed_fi, output_folder_kmean, "kmean");
	separate_bed_by_class(hclust_result, data_mat, bed_fi, output_folder_hclust, "hclust");
	separate_bed_by_class(fisher_result, data_mat, bed_fi, output_folder_fisher, "fisher");
	
}

################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);



if(length(bam_file) != 1 || length(gtf_file) != 1 || length(output_folder) != 1)
{
	stop("\n This script considers one bam file at a time. gtf file and output folder arguments should be unique.\n\n");
}

checkingOutputFolder(output_folder);

pal1 <- c("wheat1", "red3");

cat("Creating a BamFileList\n");
bamfiles <- BamFileList(bam_file);

cat("Creating the genes database from biomarRt using ensembl\n");
txdb <- makeTxDbFromGFF(gtf_file, format="gtf");

cat("Making a Granges list of all exons by genes\n");
ebg <- exonsBy(txdb, by="gene");

cat("Computing the read counts\n\n");
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
		mode="Union",
		singleEnd=single_end,
		ignore.strand= TRUE,
		fragments= if(single_end) FALSE else TRUE);


cat("Removing silent genes\n");

filtered_raw_count <- as.matrix(assay(se)[-which(assay(se) < 2),]);


cat("Plotting distribution of raw expression\n");

plot_histogram_table(output_folder, "nb_genes_raw_counts", log2(filtered_raw_count), "Raw counts");

cat("\n !!!!!!!!!!!! EdgeR process !!!!!!!!!!!!\n");

bed_fi <- read.table(bed_anno,stringsAsFactors=F)
length_vec <- bed_fi$V3 - bed_fi$V2;
names(length_vec) <- bed_fi$V4;

names_count_table <- rownames(filtered_raw_count);

indexes_for_length <- match(names_count_table, names(length_vec));

if(length(which(is.na(match(names_count_table, names(length_vec))))) != 0)
{
	stop("\n The gtf and bed files provided contain different annotations. Make sure they do.\n");
}

if(length(indexes_for_length) != length(names_count_table))
{
	stop("\n The names of annotations should be unique in the bed and gtf files\n");
}

length_vec <- length_vec[indexes_for_length];
edgeR_rpkm_matrix <- rpkm(filtered_raw_count, length_vec);
edgeR_rpkm_log_matrix <- rpkm(filtered_raw_count, length_vec, log=TRUE);

edgeR_cpm_matrix <- cpm(filtered_raw_count);
edgeR_cpm_log_matrix <- cpm(filtered_raw_count, log=TRUE);


plot_histogram_table(output_folder, "nb_genes_rpkm_counts", edgeR_rpkm_matrix, "RPKM counts");
plot_histogram_table(output_folder, "nb_genes_rpkm_counts_log", edgeR_rpkm_log_matrix, "log2(RPKM counts)");
plot_histogram_table(output_folder, "nb_genes_tpm_counts", edgeR_cpm_matrix, "TPM counts");
plot_histogram_table(output_folder, "nb_genes_tpm_counts_log", edgeR_cpm_log_matrix, "log2(TPM counts)");





##########################
# Data classification
##########################



cat("Treating intervals for rpkm matrix");

output_folder_edgR_rpkm <- paste(output_folder, "edgeR_rpkm/", sep="");
checkingOutputFolder(output_folder_edgR_rpkm);

output_folder_edgR_rpkm_log <- paste(output_folder, "edgeR_rpkm_log/", sep="");
checkingOutputFolder(output_folder_edgR_rpkm_log);

output_folder_edgR_cpm <- paste(output_folder, "edgeR_TPM/", sep="");
checkingOutputFolder(output_folder_edgR_cpm);

output_folder_edgR_cpm_log <- paste(output_folder, "edgeR_TPM_log/", sep="");
checkingOutputFolder(output_folder_edgR_cpm_log);


output_folder_raw <- paste(output_folder, "raw_counts/", sep="");
checkingOutputFolder(output_folder_raw);


computing_class_intervals(edgeR_rpkm_matrix, nb_of_groups, output_folder_edgR_rpkm, bed_fi);
computing_class_intervals(edgeR_rpkm_log_matrix, nb_of_groups, output_folder_edgR_rpkm_log, bed_fi);
computing_class_intervals(edgeR_cpm_matrix, nb_of_groups, output_folder_edgR_cpm, bed_fi);
computing_class_intervals(edgeR_cpm_log_matrix, nb_of_groups, output_folder_edgR_cpm_log, bed_fi);
computing_class_intervals(filtered_raw_count, nb_of_groups, output_folder_raw, bed_fi);

















