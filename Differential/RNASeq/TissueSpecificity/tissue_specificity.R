###############
# This script assesses the tissue specificity of a list of gene given in GFF by using different databases and the tau calculation (other metrics can be found in the tables).
# The tau metrics was developed in Yanai et al, Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification, Bioinformatics, 2005 (see methods)
# "τ values interpolate the entire range between 0 for housekeeping genes and 1 for strictly one-tissue-specific genes. It is seen (Fig. 2A) that τ values near 0 and 1 tend to be more probable than the intermediate values, generating a U-shaped distribution. However, as many as 57% of all profiles have intermediate specificities: 0.15 ≤ τ ≤ 0.85, constituting the largest group, greater than the housekeeping and one-tissue-specific sets combined."
# The tables of expression and tau metric were retrieved from the supplementary data of:
# Kryuchkova-Mostacci and Robinson-Rechavi, A benchmark of gene expression tissue-specificity metrics, Briefings in Bioinformatics, 18(2), 2017, 205–214
# In this benchmark study, it is claimed that the tau metrics is the best to use.
# Link found in the supplementary material:
#		https://figshare.com/articles/A_benchmark_of_gene_expression_tissue_specificity_metrics/1558257
# The HumBrawandTScomparisonTable_9_6h.txt is coming from:
#		Brawand et al, The evolution of gene expression levels in mammalian organs, Nature, 2011
# HumFagerbergTScomparisonTable_9_27.txt is coming from:
#		Fagerberg et al, Analysis of the Human Tissue-specific Expression by Genome-wide Integration of Transcriptomics and Antibody-based Proteomics, Mol. and cell. proteomics, 2014
# HumGeTScomparisonTable_9_32hG.txt is coming from:
#		Ge X et al. Interpreting expression profiles of cancers by genome-wide survey of breadth of expression in normal tissues. Genomics 2005;86:127–41
# The different metrics are used as follows:
# tau: housekeeping = 0, tissue specific = 1
# gini: housekeeping = 0, tissue specific = (n-1)/n
# tsi: housekeeping = 0, tissue specific = 1
# counts: housekeeping = n, tissue specific = 1
# EEi: housekeeping = 0, tissue specific = > 5
# Hg: housekeeping = log2(n), tissue specific = 0
# zscore: housekeeping = 0, tissue specific = >3
# PEM: housekeeping = 0, tissue specific = 1
# SPM: housekeeping = 0, tissue specific = 1
# Descostes, March 2017
###############


library("clusterProfiler");
library("ggplot2");
library("RColorBrewer");
library("NGSprofiling");



################
# PARAMETERS
################


expression_table_vec <- c("/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumBrawandTScomparisonTable_9_6h.txt",
		"/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumFagerbergTScomparisonTable_9_27.txt",
		"/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumGeTScomparisonTable_9_32hG.txt");

expression_table_name_vec <- c("Brawand", "Fagerberg", "Gemicroarray");

genes_files_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category1/refseq.gff",
		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category2/refseq.gff",
		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category3/refseq.gff",
		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category4/refseq.gff",
		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category5/refseq.gff",
		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category6/refseq.gff");

genes_files_name_vec <- c("category1", "category2", "category3", "category4", "category5", "category6");

species_name <- "human";

output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/tissue_specificity/heatmap_categories/"; 
col_vec <- c("#D53E4F", "#FC8D59", "#FEE08B", "#E6F598", "#99D594", "#3288BD");
output_format <- "ps";

#expression_table_vec <- c("/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumBrawandTScomparisonTable_9_6h.txt",
#		"/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumFagerbergTScomparisonTable_9_27.txt",
#		"/ifs/home/descon01/cluster/Annotations/human/hg19/tissue_expression/HumGeTScomparisonTable_9_32hG.txt");
#expression_table_name_vec <- c("Brawand", "Fagerberg", "Gemicroarray");
#genes_files_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/brd2/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/H3K27me3/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/H3K36me2/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/H3K36me3/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/hgdf2/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/ledgf/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/PolII/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/peaks_with_refseq/spt16_merge/refseq.gff");
#genes_files_name_vec <- c("brd2-genes",
#		"H3K27me3-genes",
#		"H3K36me2-genes",
#		"H3K36me3-genes",
#		"hgdf2-genes",
#		"ledgf-genes",
#		"PolII-genes",
#		"spt16_merge-genes");
#species_name <- "human";
#output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/tissue_specificity/refseq_identified_with_peaks/";





################


################
# FUNCTION
################




plot_density <- function(density_values, x_lab, comparison_title, output_folder, colors, expnames_vec, output_format)
{
	min_y <- min(unlist(lapply(density_values, function(val){return(min(val$y))})));
	max_y <- max(unlist(lapply(density_values, function(val){return(max(val$y))})));
	
	min_x <- min(unlist(lapply(density_values, function(val){return(min(val$x))})));
	max_x <- max(unlist(lapply(density_values, function(val){return(max(val$x))})));
	
	if(output_format == "png")
	{
		png(filename=paste(output_folder, comparison_title,".png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder, comparison_title, ".ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder, comparison_title,".pdf"), width=10, height=10)
 }
	
	plot(NULL,xlim= c(min_x, max_x),ylim= c(min_y, max_y), xlab= x_lab, ylab="density");
	for(i in 1:length(density_values)) 
	{
		lines(density_values[[i]], col=colors[i]);
	}
	legend("topright", expnames_vec, col=colors, lty=1)
	dev.off();
}



################


##############
# MAIN
##############


metrics_method_vec <- c("tau", "gini", "tsi", "counts", "EEi", "Hg", "zscore", "PEM", "SPM");  


if(length(expression_table_vec) != length(expression_table_name_vec))
{
	stop("A name should be given for each expression table\n");
}


if(length(genes_files_vec) != length(genes_files_name_vec))
{
	stop("A name should be given for each genes file\n");
}


#if(metrics_method != "tau" && metrics_method != "gini" && metrics_method != "tsi" && metrics_method != "counts" && metrics_method != "EEi" && metrics_method != "Hg" && metrics_method != "zscore" && metrics_method != "PEM" && metrics_method != "SPM")
#{
#	stop("\n metrics_method should be tau, gini, tsi, counts, EEi, Hg, zscore, PEM or SPM");
#}


# Reading the genes
genes_list <- lapply(genes_files_vec, function(x){return(read.table(x, stringsAsFactors=F))});


# Converting the refseq ID to ensembl

if(species_name == "human") name_database <- "org.Hs.eg.db" else name_database <-"org.Mm.eg.db";
	

different_id_list <- lapply(genes_list, function(x, database_name){
			
			result <- bitr(x$V3, fromType="REFSEQ", toType=c("ENSEMBL"), annoDb= database_name);
			return(result);
			
		}, name_database); 
		
		

for(i in 1:length(expression_table_vec)) 
{
	cat("Processing ", expression_table_name_vec[i], "\n");
	expression_data <- read.table(expression_table_vec[i], header=TRUE, stringsAsFactors=F);
	
	cat("\t Retrieving the annotation indexes\n");
	gene_indexes_list <- lapply(different_id_list, function(id_vec, expr_id){
				
				indexes <- match(id_vec$ENSEMBL, expr_id);
				
				if(length(which(is.na(indexes))) != 0)
				{
					cat("\t\t Losing ", length(which(is.na(indexes))), "/", length(indexes), " annotations\n");
					indexes <- indexes[-which(is.na(indexes))];
				}
				
				return(indexes);
				
			}, expression_data[,1]);
	
	for(j in 1:length(metrics_method_vec)) 
	{
		metrics_method <- metrics_method_vec[j];
		
		cat("\n\t With metrics method ", metrics_method, "\n");
		
		output_folder_tmp <- paste(output_folder, expression_table_name_vec[i], "/", metrics_method, "/", sep="");
		checkingOutputFolder(output_folder_tmp);
		
		cat("\t\t Retrieving the ", metrics_method, " values for each gene file\n");
		
		metrics_values_vec <- if(metrics_method == "tau")expression_data$Tau else if(metrics_method == "gini")expression_data$Gini else if(metrics_method == "tsi")expression_data$Tsi else if(metrics_method == "counts")expression_data$Counts else if(metrics_method == "EEi")expression_data$Ee else if(metrics_method == "Hg")expression_data$Hg else if(metrics_method == "zscore")expression_data$Zscore else if(metrics_method == "PEM")expression_data$Pem else if(metrics_method == "SPM")expression_data$Spm else stop("\n problem in the script regarding the metrics method\n");
			
		values_list <- lapply(gene_indexes_list, function(index_vec, values_vec){
					
					result <- values_vec[index_vec];
					return(result);
					
				}, metrics_values_vec);
		
		
		cat("\t\t Plotting the tissue specificity barplot for each category\n");
		
		hist_list <- lapply(values_list, function(x){
					
					p <- hist(x, breaks=100, plot=FALSE);
					p$counts <- p$counts/length(x);
					return(p);
				});
		y_max_vec <- c(0, max(unlist(lapply(hist_list, function(x){return(max(x$counts))}))));
		
		# All in 1
		
		if(output_format == "png")
		{
			png(filename=paste(output_folder_tmp, "all_hist.png",sep=""), width = 600, height = 600, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_tmp, "all_hist.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste0(output_folder_tmp, "all_hist.pdf"), width=10, height=10)
  }
		
		plot(hist_list[[1]], col= alpha(col_vec[1],0.5), ylim=y_max_vec, xlab="Tissue specificity index", ylab="Frequency/number of genes", main="");
		for(k in 2:length(hist_list)) 
		{
			plot(hist_list[[k]], col= alpha(col_vec[k],0.5), add=TRUE, ylim=y_max_vec);
		}
		abline(v= 0.15, col="firebrick", lty=2);
		abline(v= 0.85, col="firebrick", lty=2);
		dev.off();
		
		
		# Independently
		for(k in 1:length(hist_list)) 
		{
			if(output_format == "png")
			{
				png(filename=paste(output_folder_tmp, genes_files_name_vec[k],"-hist.png",sep=""), width = 600, height = 600, bg = "transparent");
			}else if(output_format == "ps"){
				cairo_ps(filename=paste(output_folder_tmp, genes_files_name_vec[k],"-hist.ps",sep=""), width = 7, height = 7, bg = "transparent");
			}else{
       pdf(file=paste(output_folder_tmp, genes_files_name_vec[k],"-hist.pdf",sep=""), width=10, height=10)
   }
			plot(hist_list[[k]], col= alpha(col_vec[k],0.5), ylim=y_max_vec, xlab="Tissue specificity index", ylab="Frequency/number of genes", main= genes_files_name_vec[k]);
			abline(v= 0.15, col="firebrick", lty=2);
			abline(v= 0.85, col="firebrick", lty=2);
			dev.off();
		}
		
		
		# Independently with density
		for(k in 1:length(values_list)) 
		{
			if(output_format == "png")
			{
				png(filename=paste(output_folder_tmp, genes_files_name_vec[k],"-hist-density.png",sep=""), width = 600, height = 600, bg = "transparent");
			}else if(output_format == "ps"){
				cairo_ps(filename=paste(output_folder_tmp, genes_files_name_vec[k],"-hist-density.ps",sep=""), width = 7, height = 7, bg = "transparent");
			}else{
       pdf(file=paste(output_folder_tmp, genes_files_name_vec[k],"-hist-density.pdf",sep=""), width=10, height=10)
   }
			
			hist(values_list[[k]], freq=FALSE, col= alpha(col_vec[k],0.5), xlab="Tissue specificity index", ylab="density", main= genes_files_name_vec[k])
			lines(density(values_list[[k]], from=0, to=max(values_list[[k]])), col="firebrick");
			abline(v= 0.15, col="firebrick", lty=2);
			abline(v= 0.85, col="firebrick", lty=2);
			dev.off();
		}
		
		#Plotting all density in 1
		
		density_values <- lapply(values_list, density);
		density_values_log2 <- lapply(values_list, function(x){return(density(log2(x)))});
		
		plot_density(density_values, "tissue-specificity index", "density_distribution", output_folder_tmp, col_vec, genes_files_name_vec, output_format);
		plot_density(density_values_log2, "log2(tissue-specificity index)", "density_distribution_log2", output_folder_tmp, col_vec, genes_files_name_vec, output_format);
		
		#Barplots of the number of housekeeping and tissue specific genes
		
		nb_housekeeping_norm_vec <- unlist(lapply(values_list, function(x){ return(length(which(x < 0.15))/length(x));}));
		nb_tissue_specific_norm_vec <- unlist(lapply(values_list, function(x){ return(length(which(x > 0.85))/length(x));}));
		
		if(output_format == "png")
		{
			png(filename=paste(output_folder_tmp, "housekeeping_norm.png",sep=""), width = 600, height = 600, bg = "transparent");
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_tmp, "housekeeping_norm.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(output_folder_tmp, "housekeeping_norm.pdf",sep=""), width=10, height=10)
  }
		barplot(nb_housekeeping_norm_vec, beside=TRUE, col=col_vec, names.arg= genes_files_name_vec, ylab="Nb housekeeping/nb of genes", las=2);
		dev.off();
		
		if(output_format == "png")
		{
			png(filename=paste(output_folder_tmp, "tissuespe_norm.png",sep=""), width = 600, height = 600, bg = "transparent");	
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_tmp, "tissuespe_norm.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(output_folder_tmp, "tissuespe_norm.pdf",sep=""), width=10, height=10)
  }
		barplot(nb_tissue_specific_norm_vec, beside=TRUE, col=col_vec, names.arg= genes_files_name_vec, ylab="Nb tissue specific/nb of genes", las=2);
		dev.off();
		
		nb_housekeeping_vec <- unlist(lapply(values_list, function(x){ return(length(which(x < 0.15)));}));
		nb_tissue_specific_vec <- unlist(lapply(values_list, function(x){ return(length(which(x > 0.85)));}));
		
		if(output_format == "png")
		{
			png(filename=paste(output_folder_tmp, "housekeeping.png",sep=""), width = 600, height = 600, bg = "transparent");
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_tmp, "housekeeping.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(output_folder_tmp, "housekeeping.pdf",sep=""), width=10, height=10)
  }
		barplot(nb_housekeeping_vec, beside=TRUE, col=col_vec, names.arg= genes_files_name_vec, ylab="Nb housekeeping/nb of genes", las=2);
		dev.off();
		
		if(output_format == "png")
		{
			png(filename=paste(output_folder_tmp, "tissuespe.png",sep=""), width = 600, height = 600, bg = "transparent");
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_tmp, "tissuespe.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(output_folder_tmp, "tissuespe.pdf",sep=""), width=10, height=10)
  }
		barplot(nb_tissue_specific_vec, beside=TRUE, col=col_vec, names.arg= genes_files_name_vec, ylab="Nb tissue specific/nb of genes", las=2);
		dev.off();
	}
}





