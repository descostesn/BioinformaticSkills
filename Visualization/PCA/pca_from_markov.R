###########
# This script performs a pca from objects. It is developped to perform pca on states detected with a hidden markov model.
# Descostes March 2017
###########


library("RColorBrewer");
library("NGSprofiling");
library("ggplot2");


################
# PARAMETERS
################


object_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_Flag-hdgf2_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_Flag-ledgf_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_H3K27me3_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_H3K36me2_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_H3K36me3_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_RNAPolII_spaceSizes_NA_ProLength_15000",
		"/ifs/home/descon01/analysis/fact_ledgf/objects/293/hiddenMarkovModel/model4_80kb_4st/gary_dense_bed/GenesSelect_SPT16_spaceSizes_NA_ProLength_15000");

names_object_vec <- c("Hdgf2",
		"Flag-ledgf",
		"H3K27me3",
		"H3K36me2",
		"H3K36me3",
		"RNAPolII",
		"SPT16");


original_gff_path <- c("/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/chromHMM/293TRex/after_factGary/model4_gary-window-80000/learned_states/4-states/293TRex_4_model4_gary_dense.bed");

output_folder <-"/ifs/home/descon01/analysis/fact_ledgf/hiddenMarkovModel/analysis/PCA/";

col_vec_states <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99");
		
output_format <- "pdf"
		
################




##############
# MAIN
##############


checkingOutputFolder(output_folder);



if(length(original_gff_path) != 1)
{
	stop("\n The original gff should be unique\n");
}


# Reading all objects

object_list <- lapply(object_vec, function(x){load(x); return(genesSelected)});

# Performing mean on each interval

mean_INSIDE_list <- lapply(object_list, function(x){
			values_inside <- lapply(x, "[[", "probe.valueINSIDE");
			mean_inside <- as.numeric(unlist(lapply(values_inside, mean)));
			return(mean_inside);}) 
		
matrix_inside_mean <- do.call(cbind, mean_INSIDE_list);
colnames(matrix_inside_mean) <- names_object_vec;


# Reading the gff file and computing a semicolon vec

original_gff <- read.table(original_gff_path, stringsAsFactors=F);
semicol_original_gff <- paste(original_gff$V1, original_gff$V2, original_gff$V3, sep=":");

# Checking that each object is in the same order than the original gff

result_tmp <- lapply(object_list, function(x, semicol_original){
			
			semicol_object <- paste(as.character(unlist(lapply(x,"[[","chr"))), format(unlist(lapply(x,"[[","gene.begin")), scientific=F, trim=T), format(unlist(lapply(x,"[[","gene.end")), scientific=F, trim=T), sep=":");
			index <- match(semicol_object, semicol_original);
			if(length(which(is.na(index))) != 0) stop("\n Some annotations were not retrieved\n");
			if(is.unsorted(index)) stop("\n The object and original gff are not in the same order, code ordering\n");
		}, semicol_original_gff);

# Creating the color vector corresponding to the annotations and defined by the states detected by the markov model

state_vec <- original_gff$V4;

for(i in 1:length(col_vec_states))
{
	state_vec[which(state_vec == i)] <- col_vec_states[i];
}

# Performing the pca

pca_result <- prcomp(matrix_inside_mean);

if(output_format == "png")
{
	png(filename=paste(output_folder, "variance_vs_component.png",sep=""), width = 600, height = 600, bg = "transparent")
}else if(output_format == "ps"){
	cairo_ps(filename=paste(output_folder, "variance_vs_component.ps",sep=""), width = 7, height = 7, bg = "transparent");
}else{
    pdf(file=paste(output_folder, "variance_vs_component.pdf",sep=""), width=10, height=10)
}
plot(pca_result, type = "l")
dev.off();


df_for_plot <- data.frame(pca_result$x, group= state_vec);
percentVar <- pca_result$sdev^2 / sum( pca_result$sdev^2 );
colname_vec <- colnames(df_for_plot);


cat("Plotting PCA\n");

for(i in 1:(length(colname_vec)-2))
{
	for(j in 2:(length(colname_vec)-1))
	{
		principal_componant_1 <- colname_vec[i];
		principal_componant_2 <- colname_vec[j];
		
		ggplot(data=df_for_plot, aes_string(x= principal_componant_1, y= principal_componant_2)) + geom_point(size=0.2, colour = state_vec) + 
				xlab(paste0(principal_componant_1, ": ",round(percentVar[i] * 100),"% variance")) +
				ylab(paste0(principal_componant_2, ": ",round(percentVar[j] * 100),"% variance")) +
				coord_fixed() + 
				theme( panel.background = element_rect(fill = "white", colour = "black"));
		
		if(output_format == "png")
		{
			ggsave(paste(principal_componant_1, "_", principal_componant_2, ".png",sep=""), plot = last_plot(), device = "png", path = output_folder);
		}else if(output_format == "ps"){
			ggsave(paste(principal_componant_1, "_", principal_componant_2, ".ps",sep=""), plot = last_plot(), device = "ps", path = output_folder);
		}else{
      ggsave(paste(principal_componant_1, "_", principal_componant_2, ".pdf",sep=""), plot = last_plot(), device = "pdf", path = output_folder);
  }
	}
}








