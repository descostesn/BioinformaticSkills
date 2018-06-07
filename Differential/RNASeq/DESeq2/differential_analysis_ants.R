####################
# This script aims at performing differential analysis with a non reference genome.
# Nicolas Descostes February 2016
####################


library("BiocParallel");
library("DESeq2");
library("vsn");
library("pheatmap");
library("RColorBrewer");
library("sva");
library("Rargs");
library("Rsubread");
library("ggplot2");

################
# PARAMETERS
################


# Define the required options
paramsDefinition=list();

paramsDefinition[["--bamFilesVec"]]=list(variableName="bam_files_vec", numeric=F, mandatory=T, description="A space separated list of bam files to compare.");
paramsDefinition[["--nbCPU"]]=list(variableName="nb_cpu", numeric=T, mandatory=T, description="Number of cpu to use.");
paramsDefinition[["--singleEnd"]]=list(variableName="single_end", numeric=F, mandatory=T, description="Boolean indicating if the data are single or paired ended.", postConversion=as.logical);
paramsDefinition[["--samplesInfoFile"]]=list(variableName="samples_info_file", numeric=F, mandatory=T, description="path to the file containing information about samples. Format: Sample name,cell,condition,avglength,experiment,sample,biosample");
paramsDefinition[["--nameOfReference"]]=list(variableName="name_of_reference", numeric=F, mandatory=T, description="Name of the reference condition.");
paramsDefinition[["--outputFolder"]]=list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the output folder.");
paramsDefinition[["--comparisonName"]]=list(variableName="comparison_name", numeric=F, mandatory=T, description="Name of the analysis.");
paramsDefinition[["--refseqAnno"]]=list(variableName="refseq_anno", numeric=F, mandatory=T, description="File path of the file containing annotations. This should be in gtf format");
paramsDefinition[["--expnamesVec"]]=list(variableName="expnames_vec", numeric=F, mandatory=T, description="A space separated vector of exp names corresponding to the bam files. !! Should have the same names than in sampleInfoFile !!!");

#Optional
paramsDefinition[["--singleGeneToHighlight"]]=list(variableName="single_gene_to_highlight", numeric=F, mandatory=F, description="Name of gene to highlight in the plot.", default=NA);
paramsDefinition[["--labelsSingleGene"]]=list(variableName="labels_single_gene", numeric=F, mandatory=F, description="Name of the label on the plot for the gene to highlight.", default=NA);

################



################
# FUNCTION
################

heatmapDistance <- function(object, output_folder, filename)
{
	sampleDists <- dist(t(assay(object)));
	sampleDistMatrix <- as.matrix(sampleDists);
	rownames(sampleDistMatrix) <- object$Sample.name;
	colnames(sampleDistMatrix) <- object$Sample.name;
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255);
	
	png(filename=paste(output_folder, filename,sep=""), width = 600, height = 600, bg = "transparent");
	pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors);
	dev.off();
}


checking_outputFolder <- function(output_path)
{
	if(!file.exists(output_path))
	{
		dir.create(output_path, recursive = TRUE)
	}
}


	
qualityAssessment_and_filtering <- function(DESeqDataSet, outputFolder, notAllZero)
{
	cat("\t Applying two different normalization strategies\n");
	rld <- rlog(DESeqDataSet); #Regularized log
	vsd <- varianceStabilizingTransformation(DESeqDataSet); #Variance stabilizing transformation, do not use if the differences are due to the experimental conditions
	
	cat("\t Plotting the effect of the different normalization on the variance\n");
	png(filename=paste(outputFolder, "effect_on_variance_log.png",sep=""), width = 600, height = 600, bg = "transparent");
	meanSdPlot(log2(counts(DESeqDataSet,normalized=TRUE)[notAllZero,] + 1), xlab = "log2(mean > 0 + 1)");
	dev.off();
	
	png(filename=paste(outputFolder, "effect_on_variance_rlog.png",sep=""), width = 600, height = 600, bg = "transparent");
	meanSdPlot(assay(rld[notAllZero,]), xlab = "Regularized log");
	dev.off();
	
	png(filename=paste(outputFolder, "effect_on_variance_vst.png",sep=""), width = 600, height = 600, bg = "transparent");
	meanSdPlot(assay(vsd[notAllZero,]), xlab = "Variance stabilizing transformation");
	dev.off();
	
	cat("\t Creating heatmap of the count matrix\n");
	select <- order(rowMeans(counts(DESeqDataSet,normalized=TRUE)),decreasing=TRUE)[1:20]; # normalized = divide the counts by the size factors or normalization factors before returning
	nt <- normTransform(DESeqDataSet); # defaults to log2(x+1)
	log2.norm.counts <- assay(nt)[select,];
	
	png(filename=paste(outputFolder, "heatmap_count_matrix_pseudo.png",sep=""), width = 600, height = 600, bg = "transparent");
	pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, main="log2 + 1");
	dev.off();
	
	png(filename=paste(outputFolder, "heatmap_count_matrix_rld.png",sep=""), width = 600, height = 600, bg = "transparent");
	pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, main = "rld");
	dev.off();
	
	png(filename=paste(outputFolder, "heatmap_count_matrix_vsd.png",sep=""), width = 600, height = 600, bg = "transparent");
	pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, main="vsd");
	dev.off();
	
	cat("\t Creating heatmaps of the sample to sample euclidean distance\n");
	heatmapDistance(nt, outputFolder, "heatmap_distance_pseudocounts.png");
	heatmapDistance(rld, outputFolder, "heatmap_distance_rld.png");
	heatmapDistance(vsd, outputFolder, "heatmap_distance_vsd.png");
	
	cat("\t Plotting the PCA of the sample distance");
	png(filename=paste(outputFolder, "PCA_pseudo_deseq.png",sep=""), width = 600, height = 600, bg = "transparent");
	tmp <- plotPCA(nt);
	print(tmp);
	dev.off();
	rm(tmp);
	png(filename=paste(outputFolder, "PCA_rld_deseq.png",sep=""), width = 600, height = 600, bg = "transparent");
	tmp <- plotPCA(rld);
	print(tmp);
	dev.off();
	rm(tmp);
	png(filename=paste(outputFolder, "PCA_vsd_deseq.png",sep=""), width = 600, height = 600, bg = "transparent");
	tmp <- plotPCA(vsd);
	print(tmp);
	dev.off();
	rm(tmp);
	
	# Plot PCA with ggplot

	d = plotPCA(rld, returnData=TRUE);
	ggplot(d, aes(x=PC1,y=PC2,col=condition,label=name)) + geom_point() + geom_text(aes(y = PC2 + 0.5), position = position_dodge(0.9), vjust = 0.3, hjust=0.7, angle=45);
	ggsave("PCA_rld_ggplot.png", plot = last_plot(), device = "png", path = outputFolder);
	rm(d);
	
	d = plotPCA(nt, returnData=TRUE);
	ggplot(d, aes(x=PC1,y=PC2,col=condition,label=name)) + geom_point() + geom_text(aes(y = PC2 + 0.5), position = position_dodge(0.9), vjust = 0.3, hjust=0.7, angle=45);
	ggsave("PCA_pseudo_ggplot.png", plot = last_plot(), device = "png", path = outputFolder);
	rm(d);
	
	d = plotPCA(vsd, returnData=TRUE);
	ggplot(d, aes(x=PC1,y=PC2,col=condition,label=name)) + geom_point() + geom_text(aes(y = PC2 + 0.5), position = position_dodge(0.9), vjust = 0.3, hjust=0.7, angle=45);
	ggsave("PCA_vsd_ggplot.png", plot = last_plot(), device = "png", path = outputFolder);
	rm(d);
}



barplot_and_writeTables <- function(outputFolder, title_barplot, upregulated, downregulated, DESeqResults, title_firstable, title_secondtable, title_thirdtable)
{
	png(filename=paste(outputFolder, title_barplot,sep=""), width = 600, height = 600, bg = "transparent");
	barplot(c(nrow(upregulated),nrow(downregulated)), names.arg=c(paste("upregulated: ", nrow(upregulated), sep=""), paste("downregulated: ", nrow(downregulated),sep="")));
	dev.off();
	
	cat("\t Writing the tables of the differential analysis to ", outputFolder, "\n");
	write.table(DESeqResults, file= paste(outputFolder, title_firstable, sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
	write.table(downregulated, file= paste(outputFolder, title_secondtable, sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
	write.table(upregulated, file= paste(outputFolder, title_thirdtable, sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
}




writing_differential_tables <- function(DESeqDataSet, outputFolder)
{
	# Filtering rows having a basemean at 0 and then, having a p-value equal to NA
	result_diff <- results(DESeqDataSet);
	result_diff <- result_diff[which(result_diff[,1] != 0),];
	result_diff <- result_diff[which(!is.na(result_diff$padj)),];
	
	write.table(result_diff, file=paste(outputFolder, "all-genes-diff.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
	
	cat("\t Ordering by p-values\n");
	result_diff_Ordered <- result_diff[order(result_diff$padj),];
	
	cat("\t Here is a summary of the differential analysis\n");
	summary(result_diff_Ordered, alpha=0.05);
	
	number_pvalues_inf <- sum(result_diff_Ordered$padj <= 0.05, na.rm=TRUE);
	number_total_row <- nrow(result_diff_Ordered);
	number_up_twoFold_pval <- sum(result_diff$log2FoldChange >= 1 & result_diff$padj <= 0.05, na.rm=TRUE);
	number_down_twoFold_pval <- sum(result_diff$log2FoldChange <= -1 & result_diff$padj <= 0.05, na.rm=TRUE);
	
	cat("\t Number of p-values < 0.05: ", number_pvalues_inf, "/", number_total_row, " (", round((number_pvalues_inf*100)/number_total_row), " %) of the rows\n");
	cat("\t Number of p-values < 0.05 and upregulated by at least 2 folds: ", number_up_twoFold_pval, "/", number_total_row, " (", round((number_up_twoFold_pval*100)/number_total_row), " %) of the rows\n");
	cat("\t Number of p-values < 0.05 and downregulated by at least 2 folds: ", number_down_twoFold_pval, "/", number_total_row, " (", round((number_down_twoFold_pval*100)/number_total_row), " %) of the rows\n");
	
	cat("\t Keeping only rows having a p-value < 0.05\n");
	result_diff_Ordered_pvalfiltered <- result_diff_Ordered[which(result_diff_Ordered$padj <= 0.05),];
	
	up_regulated_nonZerofold <- result_diff_Ordered_pvalfiltered[which(result_diff_Ordered_pvalfiltered$log2FoldChange > 0), ];
	down_regulated_nonZerofold <- result_diff_Ordered_pvalfiltered[which(result_diff_Ordered_pvalfiltered$log2FoldChange < 0), ];
	
	up_regulated_2fold <- result_diff_Ordered_pvalfiltered[which(result_diff_Ordered_pvalfiltered$log2FoldChange >= 1), ];
	down_regulated_2fold <- result_diff_Ordered_pvalfiltered[which(result_diff_Ordered_pvalfiltered$log2FoldChange <= -1), ];
	
	barplot_and_writeTables(outputFolder, "barplot_diff_nonzero.png", up_regulated_nonZerofold, down_regulated_nonZerofold, result_diff_Ordered_pvalfiltered, "down_and_up_regulated_significantPval_nonzerofold.txt", "downregulated_significantPval_nonzerofold.txt", "upregulated_significantPval_nonzero.txt");
	barplot_and_writeTables(outputFolder, "barplot_diff_2fold.png", up_regulated_2fold, down_regulated_2fold, result_diff_Ordered_pvalfiltered, "down_and_up_regulated_significantPval_2fold.txt", "downregulated_significantPval_2fold.txt", "upregulated_significantPval_2fold.txt");
	
	###########
	########### 
	# Keeping only 2fold with significant p-val to write the table of counts and performing a hierarchical clustering (with individual replicates, with mean and with median)
	###########
	###########
	
	cat("Writing rlog transformed counts and performing hierarchical clustering\n");
	rld <- rlog(DESeqDataSet); #Regularized log, this normalization is important to do before clustering anything
	write.table(assay(rld),file= paste(outputFolder, "rlog_counts_table_allgenes.txt", sep=""),quote=F,sep="\t");
	rlog_counts_matrix <- assay(rld);
	
	signif_diff <- c(rownames(down_regulated_2fold),rownames(up_regulated_2fold));
	index_to_keep <- match(signif_diff, rownames(rlog_counts_matrix));
	
	if(length(which(is.na(index_to_keep))) != 0)
	{
		stop("\n Problem in retrieving genes for making the hierarchical clustering\n\n");
	}
	
	signif_rlog_matrix <- assay(rld)[index_to_keep,];
	write.table(signif_rlog_matrix,file= paste(outputFolder, "signifPval_DownUP2fold_rlogCounts.txt", sep=""),quote=F,sep="\t");
	
	
	if(!is.vector(signif_rlog_matrix) && (is.matrix(signif_rlog_matrix) && nrow(signif_rlog_matrix) > 3))  #If there are some differentially expressed genes, clustering will be performed on it
	{
		#Separating the control vs ko columns
		conditions_vec <- substr(colnames(signif_rlog_matrix),1,1);
		index_WT_KO <- split(1:ncol(signif_rlog_matrix),conditions_vec);
		signif_rlog_matrix_WT <- signif_rlog_matrix[,index_WT_KO[[1]]];
		signif_rlog_matrix_KO <- signif_rlog_matrix[,index_WT_KO[[2]]];
		
		#####
		# PART1: Considering each replicate
		#####
		
		cat("\t Performing hierarchical clustering with each replicate\n");
		
		signif_rlog_matrix_norm <- signif_rlog_matrix - rowMeans(signif_rlog_matrix)
		df <- as.data.frame(colData(rld)[,c("cell","condition")]);
		outputFolder_clustering <- paste(outputFolder, "hierarchical_clustering/replicates/",sep="");
		checking_outputFolder(outputFolder_clustering);
		file_name_part1_base <- "hclust_replicates_sigPval2fold"
		
		for(nb_of_group_aggregate in c(seq(3, min(nrow(signif_rlog_matrix)-1,10),by=2),NA))
		{
			for(method_clustering in c("single", "complete", "average"))
			{
				file_name_part1 <- paste(file_name_part1_base, "_", nb_of_group_aggregate, "group_", method_clustering, sep="");
				
				png(filename=paste(outputFolder_clustering, file_name_part1, ".png",sep=""), width = 600, height = 600, bg = "transparent");
				result_clustering <- pheatmap(signif_rlog_matrix_norm, annotation_col=df, kmeans_k= nb_of_group_aggregate, clustering_method=method_clustering, treeheight_row=100, show_rownames= if(is.na(nb_of_group_aggregate)) FALSE else TRUE, main=paste(nb_of_group_aggregate, "group_", method_clustering,sep=""));
				dev.off();
				
				if(!is.na(nb_of_group_aggregate)) #Retrieving the list of annotations per group
				{
					group_list <- split(rownames(signif_rlog_matrix_norm),result_clustering[[3]]$cluster);
					catch_result <- mapply(function(groups, names, outputfolder_arg){
								
								write(groups, file=paste(outputfolder_arg, names, "list.txt", sep=""), ncolumns=1);
								
							}, group_list, as.list(paste(file_name_part1,"_group", 1:nb_of_group_aggregate,sep="")), MoreArgs=list(outputFolder_clustering))
				}
			}
		}
		
		
		#####
		# PART2: Mean values
		#####
		
		cat("\t Performing hierarchical clustering on mean\n");
		signif_rlog_matrix_WT_mean <- rowMeans(signif_rlog_matrix_WT);
		signif_rlog_matrix_KO_mean <- rowMeans(signif_rlog_matrix_KO);
		signif_rlog_matrix_WT_KO_mean <- cbind(signif_rlog_matrix_WT_mean, signif_rlog_matrix_KO_mean);
		signif_rlog_matrix_WT_KO_mean_norm <- signif_rlog_matrix_WT_KO_mean - rowMeans(signif_rlog_matrix_WT_KO_mean);
		df_mean <- df[c(1,2),];
		rownames(df_mean) <- colnames(signif_rlog_matrix_WT_KO_mean_norm);
		
		outputFolder_clustering <- paste(outputFolder, "hierarchical_clustering/by_row_mean/",sep="");
		checking_outputFolder(outputFolder_clustering);
		file_name_part2_base <- "hclust_mean_sigPval2fold"
		
		for(nb_of_group_aggregate in c(seq(3, min(nrow(signif_rlog_matrix)-1,10),by=2),NA))
		{
			for(method_clustering in c("single", "complete", "average"))
			{
				file_name_part2 <- paste(file_name_part2_base, "_", nb_of_group_aggregate, "group_", method_clustering, sep="");
				
				png(filename=paste(outputFolder_clustering, file_name_part2, ".png",sep=""), width = 600, height = 600, bg = "transparent");
				result_clustering <- pheatmap(signif_rlog_matrix_WT_KO_mean_norm, annotation_col=df_mean, kmeans_k = nb_of_group_aggregate, clustering_method=method_clustering, treeheight_row=100, show_rownames= if(is.na(nb_of_group_aggregate)) FALSE else TRUE, main=paste(nb_of_group_aggregate, "group_", method_clustering,sep=""));
				dev.off();
				
				if(!is.na(nb_of_group_aggregate)) #Retrieving the list of annotations per group
				{
					group_list <- split(rownames(signif_rlog_matrix_WT_KO_mean_norm),result_clustering[[3]]$cluster);
					catch_result <- mapply(function(groups, names, outputfolder_arg){
								
								write(groups, file=paste(outputfolder_arg, names, "list.txt", sep=""), ncolumns=1);
								
							}, group_list, as.list(paste(file_name_part2,"_group", 1:nb_of_group_aggregate,sep="")), MoreArgs=list(outputFolder_clustering))
				}
			}
		}
		
		
		#####
		# PART3: Median values
		#####
		
		cat("\t Performing hierarchical clustering on median\n");
		signif_rlog_matrix_WT_median <- rowMedians(signif_rlog_matrix_WT);
		names(signif_rlog_matrix_WT_median) <- rownames(signif_rlog_matrix_WT);
		signif_rlog_matrix_KO_median <- rowMedians(signif_rlog_matrix_KO);
		names(signif_rlog_matrix_KO_median) <- rownames(signif_rlog_matrix_KO);
		signif_rlog_matrix_WT_KO_median <- cbind(signif_rlog_matrix_WT_median, signif_rlog_matrix_KO_median);
		signif_rlog_matrix_WT_KO_median_norm <- signif_rlog_matrix_WT_KO_median - rowMeans(signif_rlog_matrix_WT_KO_median);
		df_median <- df[c(1,2),];
		rownames(df_median) <- colnames(signif_rlog_matrix_WT_KO_median_norm);
		
		outputFolder_clustering <- paste(outputFolder, "hierarchical_clustering/by_row_median/",sep="");
		checking_outputFolder(outputFolder_clustering);
		file_name_part3_base <- "hclust_median_sigPval2fold"
		
		for(nb_of_group_aggregate in c(seq(3, min(nrow(signif_rlog_matrix)-1,10),by=2),NA))
		{
			for(method_clustering in c("single", "complete", "average"))
			{
				file_name_part3 <- paste(file_name_part3_base, "_", nb_of_group_aggregate, "group_", method_clustering, sep="");
				
				png(filename=paste(outputFolder_clustering, file_name_part3, ".png",sep=""), width = 600, height = 600, bg = "transparent");
				result_clustering <- pheatmap(signif_rlog_matrix_WT_KO_median_norm, annotation_col=df_median, kmeans_k= nb_of_group_aggregate, clustering_method=method_clustering, treeheight_row=100, show_rownames= if(is.na(nb_of_group_aggregate)) FALSE else TRUE, main=paste(nb_of_group_aggregate, "group_", method_clustering,sep=""));
				dev.off();
				
				if(!is.na(nb_of_group_aggregate)) #Retrieving the list of annotations per group
				{
					group_list <- split(rownames(signif_rlog_matrix_WT_KO_median_norm),result_clustering[[3]]$cluster);
					catch_result <- mapply(function(groups, names, outputfolder_arg){
								
								write(groups, file=paste(outputfolder_arg, names, "list.txt", sep=""), ncolumns=1);
								
							}, group_list, as.list(paste(file_name_part3,"_group", 1:nb_of_group_aggregate,sep="")), MoreArgs=list(outputFolder_clustering))
				}
			}
		}
		
	}
	
	return(result_diff);
}



differential_plots <- function(outputFolder, DESeqResults, comparisonName, singleGeneToHighlight, labelsSingleGene, DESeqDataSet)
{
	cat("\t Plotting a MA plot of the fold change, pval < 0.05 are highlighted in red\n");
	png(filename=paste(outputFolder, "differential_MAplot_significantPval.png",sep=""), width = 600, height = 600, bg = "transparent");
	plotMA(DESeqResults, main= comparisonName, alpha=0.05, ylim=c(-5, 5));
	if(!is.na(singleGeneToHighlight))
	{
		with(DESeqResults[which(rownames(DESeqResults) == singleGeneToHighlight), ], {points(baseMean, log2FoldChange, col="blue", cex=1, lwd=2);text(baseMean, log2FoldChange, labels= labelsSingleGene, pos=4, col="blue");});
	}
	dev.off();
	
	cat("\t Plotting a MA plot of the fold change >= 2 and pval < 0.05 (highlighted in red)\n");
	
	df <- data.frame(mean = DESeqResults$baseMean, lfc = DESeqResults$log2FoldChange, isDE = ifelse(is.na(DESeqResults$padj), FALSE, DESeqResults$padj < 0.05 & abs(DESeqResults$log2) >= 1));
	png(filename=paste(outputFolder, "differential_MAplot_significantPval_2folds.png",sep=""), width = 600, height = 600, bg = "transparent");
	plotMA(df, main=paste(comparison_name, "-2folds",sep=""), ylim=c(-5,5));
	if(!is.na(singleGeneToHighlight))
	{
		with(DESeqResults[which(rownames(DESeqResults) == singleGeneToHighlight), ], {points(baseMean, log2FoldChange, col="blue", cex=1, lwd=2);text(baseMean, log2FoldChange, labels= labelsSingleGene, pos=4, col="blue");});
	}
	dev.off();
	
	
	cat("\t To identify the dots manually, see p 11 of the tutorial\n");
	
	cat("\t Plotting the MA plot with unshrunken values\n");
	resMLE <- results(DESeqDataSet, addMLE=TRUE);
	resMLE <- resMLE[which(resMLE[,1] != 0),];
	resMLE <- resMLE[which(!is.na(resMLE$padj)),];
	
	png(filename=paste(outputFolder, "differential_MAplot_unshrunken.png",sep=""), width = 600, height = 600, bg = "transparent");
	plotMA(data.frame(resMLE$baseMean, resMLE$lfcMLE, resMLE$padj < 0.05), main=paste(comparison_name, ": unshrunken LFC", sep=""), ylim=c(-5, 5))
	dev.off();
	
	cat("\t Plotting the counts for gene with the minimum p-value as example\n");
	png(filename=paste(outputFolder, "MAplot_minpvalue.png",sep=""), width = 600, height = 600, bg = "transparent");
	plotCounts(DESeqDataSet, gene=which.min(DESeqResults$padj)[1], intgroup="condition");
	dev.off();
	
	cat("\t Plotting the relation between p-values and mean expression\n");
	png(filename=paste(outputFolder, "pval_vs_mean.png",sep=""), width = 600, height = 600, bg = "transparent");
	qs <- c(0, quantile(DESeqResults$baseMean[DESeqResults$baseMean > 0], 0:6/6))
	bins <- cut(DESeqResults$baseMean, qs)
	levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
	ratios <- tapply(DESeqResults$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
	barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
	dev.off();
}


################


##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

checking_outputFolder(output_folder);
output_log <- paste(output_folder, "log_reports.txt", sep="");
sink(output_log);


###########
# PART 1: Reading the data
###########

		
cat("Computing the input counts matrix\n");
result_counts <- featureCounts(bam_files_vec, annot.ext= refseq_anno, isGTFAnnotationFile=TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", useMetaFeatures=TRUE, allowMultiOverlap=FALSE, isPairedEnd=FALSE)	
counts_matrix <- result_counts$counts;
colnames(counts_matrix) <- expnames_vec;

cat("Keeping only genes having more than 2 reads on average\n");
mean_nb_reads <- as.numeric(unlist(apply(counts_matrix, MARGIN=1,mean)));
to_keep <- which(mean_nb_reads > 2);
percent_kept <- (length(to_keep)*100)/nrow(counts_matrix);
cat("\t Keeping ", length(to_keep), "/", nrow(counts_matrix), "(", percent_kept, "%)");

counts_matrix <- counts_matrix[to_keep,];

cat("Reading the info file\n");
samples_info <- read.table(samples_info_file, sep="\t", header=T);
rownames(samples_info) <- samples_info$Sample.name; 

#Check if colnames of counts_matrix and rownames of sample_info are equal
test_equality <- unique(unique(rownames(samples_info)) == unique(colnames(counts_matrix)));

if(length(test_equality) != 1 && !test_equality)
{
	stop("colnames of counts_matrix and rownames of sample_info should be equal\n\n");
}

if(length(which(grep(name_of_reference, samples_info$condition) != 0)) == 0)
{
	stop("\n\n The name of your reference condition is not present in the table submitted, end of the program\n\n");
}

#Defining the nb of cpu to use
register(MulticoreParam(workers = nb_cpu,stop.on.error = TRUE, progressbar = TRUE));


cat("Creation of the DESeqDataSet from the matrix of counts\n");

ddsSE <- DESeqDataSetFromMatrix(counts_matrix, samples_info, design = ~condition);

cat("Performing the differential expression analysis\n");
ddsSE <- DESeq(ddsSE);
notAllZero <- (rowSums(counts(ddsSE))>0);

cat("Correcting for batch effect and performing the diff analysis\n");
#This code is taken from http://www.bioconductor.org/help/workflows/rnaseqGene/
dat <- counts(ddsSE, normalized=TRUE);
#Correcting for batch effect and performing the diff analysis
idx <- rowMeans(dat) > 1;
dat <- dat[idx,];
mod <- model.matrix(~ condition, colData(ddsSE));
mod0 <- model.matrix(~ 1, colData(ddsSE));
svseq <- svaseq(dat, mod, mod0, n.sv=2);

png(filename=paste(output_folder, "batch_effect.png",sep=""), width = 600, height = 600, bg = "transparent");
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ ddsSE$condition,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ ddsSE$condition,vertical=TRUE,main="SV2")
abline(h=0);
dev.off();

ddssva <- ddsSE;
ddssva$SV1 <- svseq$sv[,1];
ddssva$SV2 <- svseq$sv[,2];
design(ddssva) <- ~ SV1 + SV2 + condition;
ddssva <- DESeq(ddssva);

###########
# PART 2: Performing the quality assessment and filtering
###########

cat("Performing quality assessment and filtering without batch effect correction\n");
checking_outputFolder(paste(output_folder, "no_batchCorrection/", sep=""));
qualityAssessment_and_filtering(ddsSE, paste(output_folder, "no_batchCorrection/", sep=""), notAllZero);
cat("Performing quality assessment and filtering without batch effect correction\n");
checking_outputFolder(paste(output_folder, "batchCorrection/", sep=""));
qualityAssessment_and_filtering(ddssva, paste(output_folder, "batchCorrection/", sep=""), notAllZero);


###########
# PART 3: Writing result tables of the differential analysis
###########

cat("Writing tables for no batch correction\n");	
result_diff_noBatch <- writing_differential_tables(ddsSE, paste(output_folder, "no_batchCorrection/", sep=""));
cat("Writing tables for batch correction\n");
result_diff_Batch <- writing_differential_tables(ddssva, paste(output_folder, "batchCorrection/", sep=""));


###########
# PART 4: Exploring the data
###########


cat("Performing differential plot with no batch\n");
differential_plots(paste(output_folder, "no_batchCorrection/", sep=""), result_diff_noBatch, comparison_name, single_gene_to_highlight, labels_single_gene, ddsSE);
cat("Performing differential plot with batch\n");
differential_plots(paste(output_folder, "batchCorrection/", sep=""), result_diff_Batch, comparison_name, single_gene_to_highlight, labels_single_gene, ddssva);









