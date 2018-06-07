##########################
# This script performs boxplots from rna-seq bam files using a bed file of annotations. Three normalizations are used: DESeq geometric mean, DESeq2 regularized log and RPKM.
# Each boxplot will contain only the values different than 0.
# Descostes nov 2016
##########################


library("ggplot2")
library("Rsamtools");
library("GenomicAlignments");
library("GenomicFeatures");
library("edgeR");
library("DESeq");
library("DESeq2");
library("RColorBrewer");
library("Rargs");


################
# PARAMETERS
################


#parameters defined from the command line using RIO
paramsDefinition <- list();

# Required arguments
paramsDefinition[["--bamFileVec"]] <- list(variableName="bam_file_vec", numeric=F, mandatory=T, description="A vector of file pathes to the bam files.");
paramsDefinition[["--expnamesVec"]] <- list(variableName="expnames_vec", numeric=F, mandatory=T, description="A vector of experiments name corresponding to the bam files.");
paramsDefinition[["--stringIndexForGroups"]] <- list(variableName="string_index_for_groups", numeric=F, mandatory=T, description="Single string giving the indexes of expnames to create grouping. Each index must be separated by a '-' and each group by a '_'. Ex for creating two groups with four exp in each: 1-2-3-4_1-2-5-6");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Unique string giving the path to the output folder.");
paramsDefinition[["--refseqAnnoFile"]] <- list(variableName="refseq_anno_file", numeric=F, mandatory=T, description="A single file path to a gtf file containing refseq annotations.");
paramsDefinition[["--gffFile"]] <- list(variableName="gff_file", numeric=F, mandatory=T, description="A single file path to a gff file used to filter the refseq annotations.");
paramsDefinition[["--singleEnd"]] <- list(variableName="single_end", numeric=F, mandatory=T, description="Single logical indicating if the experiments are single end or paired-end.", postConversion=as.logical);
paramsDefinition[["--filterGffUniquely"]] <- list(variableName="filter_gff_uniquely", numeric=F, mandatory=T, description="Logical indicating if only unique annotations should be kept. Typically false to conserve all differentially expressed genes.");
paramsDefinition[["--indicateMeanWithSd"]] <- list(variableName="indicate_mean_with_sd", numeric=F, mandatory=T, description="Logical indicating if mean with standard deviation should be highlighted on the violin plot.", postConversion=as.logical);
paramsDefinition[["--indicateMean"]] <- list(variableName="indicate_mean", numeric=F, mandatory=T, description="Logical indicating if mean should be highlighted on the violin plot.", postConversion=as.logical);
paramsDefinition[["--indicateMedian"]] <- list(variableName="indicate_median", numeric=F, mandatory=T, description="Logical indicating if median should be highlighted on the violin plot.", postConversion=as.logical);
paramsDefinition[["--indicateBoxplot"]] <- list(variableName="indicate_boxplot", numeric=F, mandatory=T, description="Logical indicating if boxplot should be printed on top of the violin plot.", postConversion=as.logical);
paramsDefinition[["--indicateJitter"]] <- list(variableName="indicate_jitter", numeric=F, mandatory=T, description="Logical indicating if dots should be drawn on top of the violin plot. The position of dots will be randomized to enhancer visualization", postConversion=as.logical);
paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Numeric indicating the nb of CPU to use");
paramsDefinition[["--outputFormat"]] <- list(variableName="output_format", numeric=F, mandatory=T, description="output format which can be png, pdf or ps.");




#bam_file_vec <- c("/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/WT-D0-rep1_SORTED_PICARD_COOR.bam",
#                "/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/WT-D0-rep2_SORTED_PICARD_COOR.bam",
#                "/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/WT-D6-rep1_SORTED_PICARD_COOR.bam",
#                "/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/WT-D6-rep2_SORTED_PICARD_COOR.bam",
#                "/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/HGDF2_KO_2-D6-rep1_SORTED_PICARD_COOR.bam",
#                "/ifs/home/descon01/data/data_may2018/bam_files/gary_factLedgf_rnaseq/HGDF2_KO_2-D6-rep2_SORTED_PICARD_COOR.bam");
#expnames_vec <- c("WT-D0-rep1",
#"WT-D0-rep2",
#"WT-D6-rep1",
#"WT-D6-rep2",
#"HGDF2_KO_2-D6-rep1",
#"HGDF2_KO_2-D6-rep2");
#string_index_for_groups <- "1-1-2-2-3-3"
#output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/boxplots/rnaseqmay2018_comparison/based_upregulatedgenes_from_Day0ToDay6_downdKODay6/DESEQ2/") 
#refseq_anno_file <- "/ifs/home/descon01/cluster/Annotations/human/hg19/gtf_files/refseq_hg19.gtf";
#gff_file <- c("/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/genes_rnaseq_may2018/DESEQ2/WTD0D6_vs_WTD6KO2D6/peaks_per_circle/WTD0D6up_WTD6KO2D6down_original.gff");
#single_end <- TRUE;
#filter_gff_uniquely <- TRUE;
#indicate_mean_with_sd <- FALSE;
#indicate_mean <- TRUE;
#indicate_median <- TRUE;
#indicate_boxplot <- FALSE;
#indicate_jitter <- FALSE;
#nb_cpu <- 1;
#output_format <- "png"

################



################
# FUNCTION
################


perform_ggplot <- function(data_frame_input, norm_name, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, output_folder, outliers, output_format)
{
	g <- ggplot(data_frame_input, aes(expnames, values, color=factor(expnames))) + geom_violin(trim=FALSE) + ggtitle('Expression') + coord_flip() + theme_classic();
	g <- g + theme(panel.background=element_rect(fill='white', colour='black'),plot.title=element_text(size=20,face="bold",margin=margin(10,0,10,0)), axis.text.x=element_text(angle=90, size=10));
	#g <- g + labs(x="experiment", y=expression(paste(norm_name, " normalized expression", sep="")));
	g <- g + labs(x="experiment", y= paste(norm_name, " normalized expression", sep=""));
	g <- g + scale_color_brewer(palette="Set2", name="Experiments");
	
	if(indicate_mean_with_sd)
	{
		g <- g + stat_summary(fun.data = mean_sdl, geom="pointrange", shape=23, size=2, color="black");
	}
	
	if(indicate_mean)
	{
		g <- g + stat_summary(fun.y = mean, geom="point", shape=23, size=2, color="black");
	}
	
	if(indicate_median)
	{
		g <- g + stat_summary(fun.y = median, geom="point", shape=23, size=2);
	}
	
	
	if(indicate_boxplot)
	{
		g <- g + geom_boxplot(width=0.1, outlier.colour="NA");
	}
	
	
	if(indicate_jitter)
	{
		g <- g + geom_jitter(alpha=0.5, position= position_jitter(width = 0.1));
	}
	
 output_file_name <- paste("violinplot", norm_name, if(indicate_mean_with_sd) "meanSD", if(indicate_mean) "mean", if(indicate_median) "median", if(indicate_boxplot) "boxplot", if(indicate_jitter) "jitter", if(!outliers) "noOutliers", sep="_");
 
 if(output_format == "png")
 {
     ggsave(paste(output_file_name, ".png", sep=""), plot = last_plot(), device = "png", path = output_folder);
 }else if(output_format == "ps"){
     
     ggsave(paste(output_file_name, ".ps", sep=""), plot = last_plot(), device = "ps", path = output_folder);
 }else{
     ggsave(paste(output_file_name, ".pdf", sep=""), plot = last_plot(), device = "pdf", path = output_folder);
 }
 
}



perform_boxplot_and_ggplot <- function(expnames_for_grouping_list, input_matrix, output_folder, name_norm, colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format)
{
	
	cat("\t Formatting the input matrix\n");
	
	#Creating the matrix according to the names of experiments
	
	indexes_col_newMat <- unlist(lapply(expnames_for_grouping_list, function(x, refnames){ return(match(x,refnames))}, colnames(input_matrix)));
	
	new_input_matrix <- input_matrix[,indexes_col_newMat];
	
	
	#boxplot
	
	cat("\t Performing boxplots\n");
	
	if(output_format == "png")
 {
     png(filename=paste(output_folder, name_norm, "-boxplots_withNotches.png",sep=""), width = 600, height = 600, bg = "transparent")
 }else if(output_format == "ps"){
     cairo_ps(filename=paste(output_folder, name_norm, "-boxplots_withNotches.ps",sep=""), width = 7, height = 7, bg = "transparent");
 }else{
     pdf(file=paste(output_folder, name_norm, "-boxplots_withNotches.pdf",sep=""), width=10, height=10)
 }
 par(mfrow=c(1,2))
	boxplot(new_input_matrix, ylab= paste(name_norm, " normalized expression", sep=""), main= name_norm, col= colors_tab, outline=FALSE, at = locationBoxplots, names = unlist(expnames_for_grouping_list), las=2, cex.axis=0.7, notch=TRUE);
	boxplot(new_input_matrix, ylab= paste(name_norm, " normalized expression", sep=""), main= name_norm, col= colors_tab, outline=TRUE, at = locationBoxplots, names = unlist(expnames_for_grouping_list), las=2, cex.axis=0.7, notch=TRUE);
	dev.off();
	
 
	
 if(output_format == "png")
 {
     png(filename=paste(output_folder, name_norm, "-boxplots_noNotches.png",sep=""), width = 600, height = 600, bg = "transparent")	
 }else if(output_format == "ps"){
     cairo_ps(filename=paste(output_folder, name_norm, "-boxplots_noNotches.ps",sep=""), width = 7, height = 7, bg = "transparent");
 }else{
     pdf(file=paste(output_folder, name_norm, "-boxplots_noNotches.pdf",sep=""), width=10, height=10)
 }
	par(mfrow=c(1,2))
	boxplot(new_input_matrix, ylab= paste(name_norm, " normalized expression", sep=""), main= name_norm, col= colors_tab, outline=FALSE, at = locationBoxplots, names = unlist(expnames_for_grouping_list), las=2, cex.axis=0.7, notch=FALSE);
	boxplot(new_input_matrix, ylab= paste(name_norm, " normalized expression", sep=""), main= name_norm, col= colors_tab, outline=TRUE, at = locationBoxplots, names = unlist(expnames_for_grouping_list), las=2, cex.axis=0.7, notch=FALSE);
	dev.off();
	
	cat(" Writing the matrix\n");
	write.table(new_input_matrix, file=paste(output_folder, name_norm, "-matrix.txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE);
	
	##### Making the violin plot
	
	cat("\t Creating data frames for violin plots\n");
	names_col_mat <- rep(colnames(new_input_matrix), rep(nrow(new_input_matrix), ncol(new_input_matrix)));
	vector_fromMatrix <- as.vector(new_input_matrix);
	
	df_fromMatrix <- data.frame(expnames = names_col_mat, values = vector_fromMatrix);
	
	
	#removing outliers
	
	outliers_vec <- unlist(apply(new_input_matrix, MARGIN=2, function(x){ result <- boxplot(x, plot=FALSE); return(result$out); }))
	index_to_remove <- unique(unlist(mcmapply(function(x, values){return(which(values == x))}, outliers_vec, MoreArgs = list(df_fromMatrix$values), mc.cores = nb_cpu)));
	
	if(length(which(is.na(index_to_remove))) != 0)
	{
		stop("\n When removing outliers, NA values were found, pb in the script\n\n");
	}
	
 if(!is.null(index_to_remove))
 {
     df_fromMatrix_nooutliers <- df_fromMatrix[-index_to_remove,];
 }
	else{
     df_fromMatrix_nooutliers <- df_fromMatrix;
 }
	cat("\t Performing violin plots\n");
	
	perform_ggplot(df_fromMatrix_nooutliers, name_norm, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, output_folder, FALSE, output_format);
	perform_ggplot(df_fromMatrix, name_norm, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, output_folder, TRUE, output_format);
	
}


################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


expnames_for_grouping_list <- lapply(strsplit(strsplit(string_index_for_groups, "_")[[1]],"-"), function(x){return(expnames_vec[as.numeric(x)])});


if(!file.exists(output_folder))
{
	dir.create(output_folder, recursive = TRUE)
}


## Checking parameters


if(length(bam_file_vec) != length(expnames_vec))
{
	stop("\n bam files vec and exp names vec do not have the same length\n");
}

if(length(output_folder) != 1)
{
	stop("\n The path to the output folder should be unique.\n\n");
}

if(length(refseq_anno_file) != 1 || strsplit(basename(refseq_anno_file), "\\.")[[1]][2] != "gtf")
{
	stop("\n The refseq file should be unique and in gtf format.\n\n");
}

if(length(gff_file) != 1 || strsplit(basename(gff_file), "\\.")[[1]][2] != "gff")
{
	stop("\n The gff file should be unique and in gff format.\n\n");
}

if(output_format != "png" && output_format != "ps" && output_format != "pdf")
{
    stop("output_format should be ps, pdf or png\n");
}


### Starting to retrieve genes counts

cat("Creating a BamFileList\n");
bamfiles <- BamFileList(bam_file_vec);

cat("Creating the genes database from biomarRt using ensembl\n");
txdb <- makeTxDbFromGFF(refseq_anno_file, format="gtf");


cat("Making a Granges list of all exons by genes\n");
ebg <- exonsBy(txdb, by="gene");

cat("Computing the read counts\n\n");
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
		mode="Union",
		singleEnd=single_end,
		ignore.strand= TRUE,
		fragments= if(single_end) FALSE else TRUE);

colnames(se) <- expnames_vec;

raw_counts_table <- assay(se);


### Filtering the raw counts with the input gff file

cat("Filtering the gtf file\n");

gff_anno <- read.table(gff_file, stringsAsFactors=F);

if(filter_gff_uniquely)
{
	if(length(which(duplicated(gff_anno$V3))) != 0)
	{
		gff_anno <- gff_anno[-which(duplicated(gff_anno$V3)), ]
	}
}

names_gff <- gff_anno$V3;
names_count_table <- rownames(raw_counts_table);

index_in_gff <- match(names_gff, names_count_table); 

if(length(which(is.na(index_in_gff))) != 0)
{
	index_in_gff <- index_in_gff[-which(is.na(index_in_gff))];
}

raw_counts_table <- raw_counts_table[index_in_gff,];


#Checking the reverse relationship to remove genes in the gff that were not found
names_count_table <- rownames(raw_counts_table);
index_in_gtf <- match(names_count_table, names_gff);

if(length(which(is.na(index_in_gtf))) != 0) index_in_gtf <- index_in_gtf[-which(is.na(index_in_gtf))];

gff_anno <- gff_anno[index_in_gtf,];

names_gff <- gff_anno$V3;


if(length(which(is.na(match(names_count_table, names_gff)))) != 0)
{
	stop("\n Pb in the script, count table and gff should have the same elements\n");
}


if(is.unsorted(match(names_count_table, names_gff)))
{
	stop("\n Pb in the script, the count table and the gff file are not in the same order\n");
}


length_vec <- gff_anno$V5 - gff_anno$V4;
names(length_vec) <- gff_anno$V3;



##########
# PART1: EdgeR normalization
##########


cat("\n !!!!!!!!!!!! EdgeR process !!!!!!!!!!!!\n");

edgeR_rpkm_matrix <- rpkm(raw_counts_table, length_vec);
edgeR_rpkm_matrix_log2 <- rpkm(raw_counts_table, length_vec, log=TRUE);

#Removing the rows having all values equal to zero

edgeR_sum_vector <- apply(edgeR_rpkm_matrix, MARGIN=1, sum);
edgeR_sum_vector_log2 <- apply(edgeR_rpkm_matrix_log2, MARGIN=1, sum);

if(length(which(edgeR_sum_vector == 0)) != 0)
{
	edgeR_rpkm_matrix <- edgeR_rpkm_matrix[-which(edgeR_sum_vector == 0),]
}


if(length(which(edgeR_sum_vector_log2 == 0)) != 0)
{
	edgeR_rpkm_matrix_log2 <- edgeR_rpkm_matrix_log2[-which(edgeR_sum_vector_log2 == 0),]
}


##########
# PART2: DESeq normalization
##########

cat("\n !!!!!!!!!!!! DEseq process !!!!!!!!!!!!\n");

cds <- newCountDataSet(raw_counts_table, rep("ctrl", ncol(raw_counts_table)));
cds <- estimateSizeFactors(cds);
deseq_norm_matrix <- counts(cds, normalized=TRUE );




##########
# PART3: DESeq2 rlog normalization
##########

cat("\n !!!!!!!!!!!! DEseq2 rlog process !!!!!!!!!!!!\n");

deseq2_rlog_matrix <- rlog(raw_counts_table, blind = TRUE); #From the doc, "blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made. blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis.



##########
# PART4: TPM normalization
##########

cat("### TPM normalization\n")

tpm_matrix <- cpm(raw_counts_table);
tpm_matrix_log2 <- cpm(raw_counts_table, log=TRUE);


##########
# PART4: Performing the boxplot
##########



# Creating the location vector and the color vector

if(length(expnames_for_grouping_list) > 1)
{
	nb_element_pergroup <- unlist(lapply(expnames_for_grouping_list,length));		
	matrix_coor <- sapply(nb_element_pergroup, seq_len);
	
	for(i in 2:ncol(matrix_coor))
	{
		matrix_coor[,i] <- matrix_coor[,i]+ (max(matrix_coor[,i-1])+1)
	}
	
	locationBoxplots <- as.numeric(matrix_coor);
}else{
	locationBoxplots <- 1:length(expnames_for_grouping_list[[1]]);
}

colors_tab <- brewer.pal(length(locationBoxplots),"Dark2");

cat("Processing RPKM normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, edgeR_rpkm_matrix, output_folder, "RPKM", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);

cat("Processing log2(RPKM) normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, edgeR_rpkm_matrix_log2, output_folder, "log2_RPKM", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);


cat("Processing DESeq normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, deseq_norm_matrix, output_folder, "deseq", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);

cat("Processing DESeq2 normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, deseq2_rlog_matrix, output_folder, "deseq2", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);


cat("Processing TPM normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, tpm_matrix, output_folder, "TPM", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);

cat("Processing log2(TPM) normalized values\n");
perform_boxplot_and_ggplot(expnames_for_grouping_list, tpm_matrix_log2, output_folder, "log2_TPM", colors_tab, locationBoxplots, indicate_mean_with_sd, indicate_mean, indicate_median, indicate_boxplot, indicate_jitter, nb_cpu, output_format);


