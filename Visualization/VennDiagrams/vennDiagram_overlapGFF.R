###############
# This script aimns at representing differing and common intervals for a maximum of 4 gff files.
# Descostes November 2015
# update: feb 2018
###############


library("GenomicFeatures");
library("NGSprofiling");
library("ChIPpeakAnno");
library("Rargs");
library("rtracklayer");
library("biomaRt");
library("annotate");
library("gplots");
library("VennDiagram");



#################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--gffFileVec"]] <- list(variableName="gff_file_vec", numeric=F, mandatory=T, description="Vector containing file path in gff to the exp to compare.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="single path to the output folder.");
paramsDefinition[["--comparisonTitle"]] <- list(variableName="comparison_title", numeric=F, mandatory=T, description="Single string of the title of the analysis.");
paramsDefinition[["--expnamesTab"]] <- list(variableName="expnames_tab", numeric=F, mandatory=T, description="Vector containing name of each experiment given as gff files.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="Version of the genome used to retrieve features.");
paramsDefinition[["--organismName"]] <- list(variableName="organism_name", numeric=F, mandatory=T, description="String containing the name of the organism. Possible values are mouse or human.");
paramsDefinition[["--centerCriterion"]] <- list(variableName="center_criterion", numeric=F, mandatory=T, description="The center can be either the middle of the interval (coordinates) or the maximum value (max).");
paramsDefinition[["--colVec"]] <- list(variableName="col_vec", numeric=F, mandatory=T, description="Vector of colors to use for the venn diagrams.");
paramsDefinition[["--outputFormat"]] <- list(variableName="output_format", numeric=F, mandatory=T, description="output format which can be png, pdf or ps.");

#Optional argument
paramsDefinition[["--refExpMaxCentering"]] <- list(variableName="ref_exp_max_centering", numeric=F, mandatory=F, description="The reference exp name to center on max.", default=NA);
paramsDefinition[["--upstreamHeatmap"]] <- list(variableName="upstream_heatmap", numeric=T, mandatory=F, description="Number of nucleotides upstream the center of the peak on the heatmap.", default=NA);
paramsDefinition[["--downstreamHeatmap"]] <- list(variableName="downstream_heatmap", numeric=T, mandatory=F, description="Number of nucleotides downstream the center of the peak on the heatmap.", default=NA);
paramsDefinition[["--bigwigVec"]] <- list(variableName="bigwig_vec", numeric=F, mandatory=F, description="Vector of big wig files that will be used to generate heatmaps of signal.", default=NA);
paramsDefinition[["--bigwigNameVec"]] <- list(variableName="bigwig_name_vec", numeric=F, mandatory=F, description="Vector of big wig file names.", default=NA);


#gff_file_vec <- c("/ifs/home/descon01/analysis/K27M_project/overlap/peaks/data_feb2018/EZH2_8_timecourseK27M/peaks_per_circle/EZH2_0_EZH2_8_6H_EZH2_8_12H_EZH2_8_24H_EZH2_8_72H.gff", 
#        "/ifs/home/descon01/analysis/K27M_project/overlap/peaks/data_feb2018/HA-K27M/peaks_per_circle/HA_12hr-K27M_HA_24hr-K27M_HA_6hr-K27M_HA_72H_K27M_HA_K27M0Hr.gff")
#output_folder <- c("/ifs/home/descon01/analysis/K27M_project/overlap/peaks/data_feb2018/EZH2_commonOverlapK27M_vs_HAcommonOverlapinK27M-no0H/")
#comparison_title <- c("EZH2_commonOverlapK27M_vs_HAcommonOverlapinK27M-no0H")
#expnames_tab <- c("EZH2_commonOverlapK27M", "HAcommonOverlapinK27M-no0H")
#genome_version <- c("hg19");
#organism_name <- "human";
#center_criterion <- "max";
#col_vec <- c("darkblue", "darkgreen")
#bigwig_vec <- NA
#
#
#
#
#gff_file_vec <- c("/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_12hr-H3K27M-293TREX_analysis_100000.gff", 
#        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_24hr-H3K27M-293TREX_analysis_100000.gff",
#        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_6hr-H3K27M-293TREX_analysis_100000.gff",
#        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_72H_H33K27M_analysis_100000.gff",
#        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/K27me3_K27M0Hr_analysis_100000.gff")
#output_folder <- c("/ifs/home/descon01/analysis/K27M_project/overlap/peaks/data_feb2018/K27me3_timecourseK27M/")
#comparison_title <- c("K27me3_timecourseK27M")
#expnames_tab <- c("K27me3-12H K27me3-24H K27me3-6H K27me3-72H K27me3-0H")
#genome_version <- c("hg19");
#organism_name <- "human";
#center_criterion <- "max";
#col_vec <- c(paste("darkblue", "darkgreen", "darkred", "darkorange", "darkcyan", sep=" "))
#output_format <- "pdf";




################


################
# FUNCTION
################


BH_filtering_and_sorting <- function(go_results){

	results <- lapply(go_results, function(x){
		
		index_to_remove <- which(x$BH.adjusted.p.value > 0.01);
		
		if(length(index_to_remove) != 0)
		{
			x <- x[-index_to_remove, ];
		}
		
		return(x[order(x$BH.adjusted.p.value),]);
		
	});

	return(results);
};


retrieve_peaks_number <- function(peak_list, label_name)
{
	peak_number <- lapply(peak_list, length)[[label_name]];
	
	if(is.null(peak_number)) return(0) else return(peak_number);
}

################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


if(length(gff_file_vec) > 5)
{
	stop("\n This script takes at most 5 gff files as input\n");
}

if(length(gff_file_vec) < 2)
{
	stop("\n This script takes a minimum of 2 gff files as input\n");
}

if(length(gff_file_vec) != length(expnames_tab))
{
	stop("\n The number of exp names should be equal to the number of gff files");
}

if(!is.na(bigwig_vec))
{
	if(length(bigwig_name_vec) != length(bigwig_vec))
	{
		stop("\n The number of bigwig names should be equal to the number of bigwig files\n");
	}
}

if(organism_name != "mouse" && organism_name != "human")
{
	stop("\n The organism name should be mouse or human\n");
}

if(center_criterion != "coordinates" && center_criterion != "max")
{
	stop("\n Center criterion should be 'coordinates' or 'max'\n");
}

if(!is.na(bigwig_vec))
{
	if(length(grep(ref_exp_max_centering, bigwig_name_vec)) == 0 || length(ref_exp_max_centering) > 1)
	{
		stop("\n The reference name to center on max should be one of the experiment name AND bigwig name\n");
	}
}

if(!is.na(bigwig_vec) && center_criterion == "max" && is.na(ref_exp_max_centering))
{
	stop("\n You should use a reference experiment to center on max. This name has to be present in expnames AND bigwig_names\n");
}

if(organism_name == "mouse")
{
	library("org.Mm.eg.db");
}else{
	library("org.Hs.eg.db");
}


if(output_format != "ps" && output_format != "png" && output_format != "pdf")
{
	stop("output_format should be png, pdf or ps\n");
}

checkingOutputFolder(output_folder);


cat("Reading GFF files and converting to GRanges\n");

gff_list <- list();
gff_GRanges_list <- list();

for(i in 1:length(gff_file_vec)) 
{
	
	current_gff <- read.table(gff_file_vec[i], stringsAsFactors=FALSE);
		
	semicolon_feature <- paste(current_gff[,1], current_gff[,4], current_gff[,5], sep=";");
	index_duplicated <- which(duplicated(semicolon_feature));
	
	if(length(index_duplicated) != 0)
	{
		cat("\t Removing duplicated features: ", (length(index_duplicated)*100)/nrow(current_gff), "\n");
		current_gff <- current_gff[-index_duplicated,];
	}
	gff_list[[i]] <- current_gff; 
	gff_GRanges_list[[i]] <- GRanges(seqnames= gff_list[[i]][,1], ranges=IRanges(start= gff_list[[i]][,4], end = gff_list[[i]][,5], names = gff_list[[i]][,2]), strand=gff_list[[i]][,7])
}


cat("Performing the overlap\n")

if(length(gff_GRanges_list) == 2)
{
	peaks1 <- gff_GRanges_list[[1]];
	peaks2 <- gff_GRanges_list[[2]];
	
	ol <- findOverlapsOfPeaks(peaks1, peaks2, maxgap=1, connectedPeaks="keepAll");
	
}else if(length(gff_GRanges_list) == 3)
{
	peaks1 <- gff_GRanges_list[[1]];
	peaks2 <- gff_GRanges_list[[2]];
	peaks3 <- gff_GRanges_list[[3]];
	
	ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, maxgap=1, connectedPeaks="keepAll");
		
}else if(length(gff_GRanges_list) == 4)
{
	peaks1 <- gff_GRanges_list[[1]];
	peaks2 <- gff_GRanges_list[[2]];
	peaks3 <- gff_GRanges_list[[3]];
	peaks4 <- gff_GRanges_list[[4]];
	
	ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, peaks4, maxgap=1, connectedPeaks="keepAll");
}else if(length(gff_GRanges_list) == 5)
{
    peaks1 <- gff_GRanges_list[[1]];
    peaks2 <- gff_GRanges_list[[2]];
    peaks3 <- gff_GRanges_list[[3]];
    peaks4 <- gff_GRanges_list[[4]];
    peaks5 <- gff_GRanges_list[[5]];
    
    ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, peaks4, peaks5, maxgap=1, connectedPeaks="keepAll");
}

cat("Making the venn diagram");

nb_test <- max(unlist(lapply(gff_list, nrow))) + max(unlist(lapply(gff_list, nrow)));

if(output_format == "png"){
	png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno.png",sep=""), width = 600, height = 600, bg = "transparent")
}else if(output_format == "ps"){
	cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno.ps",sep=""), width = 7, height = 7, bg = "transparent");
}else{
    pdf(file=paste0(output_folder, comparison_title,"overlap-chippeakanno.pdf"), width=10, height=10)
}
result <- makeVennDiagram(ol, NameOfPeaks= expnames_tab, totalTest=nb_test, ignore.strand=TRUE, connectedPeaks="keepAll");
dev.off();

write.table(result$p.value, file=paste0(output_folder, "pvalue_stats.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE);

# Making the venn diagram venneuler


if(length(gff_GRanges_list) == 2)
{
	
	area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
	area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
	area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
	
	if(output_format == "png"){
		png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder, comparison_title,"overlap-chippeakanno-formatted.pdf"), width=10, height=10)
 }
	draw.pairwise.venn(area1 = area1 + area1_area2, 
			area2 = area2 + area1_area2, 
			cross.area = area1_area2, 
			category = expnames_tab,
			euler.d = TRUE, 
			scaled = TRUE, 
			ext.text = FALSE, 
			col = col_vec, 
			fill = col_vec, 
			cex = rep(2, 3), 
			cat.cex = rep(2.5, 2));
	dev.off();
	
}else if(length(gff_GRanges_list) == 3)
{
	
	area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
	area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
	area3 <- retrieve_peaks_number(ol$peaklist, "peaks3");
	area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
	area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3");
	area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3");
	area1_area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3");
	
	if(output_format == "png"){
		png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder, comparison_title,"overlap-chippeakanno-formatted.pdf"), width=10, height=10)
 }
	draw.triple.venn(area1 = area1 + area1_area2 + area1_area3 + area1_area2_area3, 
			area2 = area2 + area1_area2 + area2_area3 + area1_area2_area3, 
			area3 = area3 + area1_area3 + area2_area3 + area1_area2_area3, 
			n12 = area1_area2 + area1_area2_area3, 
			n23 = area2_area3 + area1_area2_area3, 
			n13 = area1_area3 + area1_area2_area3, 
			n123 = area1_area2_area3,
			category = expnames_tab,
			euler.d = TRUE, 
			scaled = TRUE, 
			ext.text = FALSE,
			col = col_vec, 
			fill = col_vec,
			cex = rep(2, 7), 
			cat.cex = rep(2.5, 3));
	dev.off();
	
}else if(length(gff_GRanges_list) == 4){
	
	area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
	area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
	area3 <- retrieve_peaks_number(ol$peaklist, "peaks3");
	area4 <- retrieve_peaks_number(ol$peaklist, "peaks4");
	
	area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
	area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3");
	area1_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4");
	area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3");
	area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4");
	area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4");
	
	area1_area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3");
	area1_area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4");
	area1_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4");
	area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4");
	
	area1_area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4");
	
	if(output_format == "png"){
		png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder, comparison_title,".pdf"), width=10, height=10)
 }
	draw.quad.venn(area1 = area1 + area1_area2 + area1_area3 + area1_area4 + area1_area2_area3 + area1_area2_area4 + area1_area3_area4 + area1_area2_area3_area4, 
			area2 = area2 + area1_area2 + area2_area3 + area2_area4 + area1_area2_area3 + area1_area2_area4 + area2_area3_area4 + area1_area2_area3_area4, 
			area3 = area3 + area1_area3 + area2_area3 + area3_area4 + area1_area2_area3 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
			area4 = area4 + area1_area4 + area2_area4 + area3_area4 + area1_area2_area4 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
			n12 = area1_area2 + area1_area2_area3 + area1_area2_area4 + area1_area2_area3_area4, 
			n13 = area1_area3 + area1_area3_area4 + area1_area2_area3 + area1_area2_area3_area4, 
			n14 = area1_area4 + area1_area2_area4 + area1_area3_area4 + area1_area2_area3_area4 , 
			n23 = area2_area3 + area1_area2_area3 + area2_area3_area4 + area1_area2_area3_area4, 
			n24 = area2_area4 + area1_area2_area4 + area2_area3_area4 + area1_area2_area3_area4,
			n34 = area3_area4 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
			n123 = area1_area2_area3 + area1_area2_area3_area4, 
			n124 = area1_area2_area4 + area1_area2_area3_area4, 
			n134 = area1_area3_area4 + area1_area2_area3_area4, 
			n234 = area2_area3_area4 + area1_area2_area3_area4, 
			n1234 = area1_area2_area3_area4, 
			category = expnames_tab, 
			euler.d = TRUE, 
			scaled = TRUE, 
			ext.text = FALSE,
			col = col_vec, 
			fill = col_vec,
			cex = rep(2, 15), 
			cat.cex = rep(2.5, 4));
	dev.off();
}else{
    
    
    
    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3");
    area4 <- retrieve_peaks_number(ol$peaklist, "peaks4");
    area5 <- retrieve_peaks_number(ol$peaklist, "peaks5");
    
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3");
    area1_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4");
    area1_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks5");
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3");
    area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4");
    area2_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks5");
    area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4");
    area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks5");
    area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks4///peaks5");
    
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3");
    area1_area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4");
    area1_area2_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks5");
    area1_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4");
    area1_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks5");
    area1_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4///peaks5");
    area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4");
    area2_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks5");
    area2_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4///peaks5");
    area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4///peaks5");
        
    area1_area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4");
    area1_area2_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks5");
    area1_area2_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4///peaks5");
    area1_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4///peaks5");
    area2_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4///peaks5");
    
    area1_area2_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4///peaks5");
    
    if(output_format == "png"){
        png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 600, height = 600, bg = "transparent")
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste0(output_folder, comparison_title,".pdf"), width=10, height=10)
    }
    
    draw.quintuple.venn(
            area1 = area1 + area1_area2 + area1_area3 + area1_area4 + area1_area5 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area1_area3_area4 + area1_area3_area5 + area1_area4_area5 + area1_area2_area3_area4 + 
                    area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area2 = area2 + area1_area2 + area2_area3 + area2_area4 + area2_area5 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area2_area3_area4 + area2_area3_area5 + area2_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area3_area5 + area1_area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area3 = area3 + area1_area3 + area2_area3 + area3_area4 + area3_area5 + area1_area2_area3 + area1_area3_area4 + area1_area3_area5 + area2_area3_area4 + area2_area3_area5 + area3_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area3_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area4 = area4 + area1_area4 + area2_area4 + area3_area4 + area4_area5 + area1_area2_area4 + area1_area3_area4 + area1_area4_area5 + area2_area3_area4 + area2_area4_area5 + area3_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            area5 = area5 + area1_area5 + area2_area5 + area3_area5 + area4_area5 + area1_area2_area5 + area1_area3_area5 + area1_area4_area5 + area2_area3_area5 + area2_area4_area5 + area3_area4_area5 + area1_area2_area3_area5 +
                    area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n12 = area1_area2 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5, 
            n13 = area1_area3 + area1_area2_area3 + area1_area3_area4 + area1_area3_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n14 = area1_area4 + area1_area2_area4 + area1_area3_area4 + area1_area4_area5 + area1_area2_area3_area4 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n15 = area1_area5 + area1_area2_area5 + area1_area3_area5 + area1_area4_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n23 = area2_area3 + area1_area2_area3 + area2_area3_area4 + area2_area3_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n24 = area2_area4 + area1_area2_area4 + area2_area3_area4 + area2_area4_area5 + area1_area2_area3_area4 + area1_area2_area4_area5  + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n25 = area2_area5 + area1_area2_area5 + area2_area3_area5 + area2_area4_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n34 = area3_area4 + area1_area3_area4 + area2_area3_area4 + area3_area4_area5 + area1_area2_area3_area4 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n35 = area3_area5 + area1_area3_area5 + area2_area3_area5 + area3_area4_area5 + area1_area2_area3_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n45 = area4_area5 + area1_area4_area5 + area2_area4_area5 + area3_area4_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n123 = area1_area2_area3 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area2_area3_area4_area5, 
            n124 = area1_area2_area4 + area1_area2_area3_area4 + area1_area2_area4_area5 + area1_area2_area3_area4_area5, 
            n125 = area1_area2_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n134 = area1_area3_area4 + area1_area2_area3_area4 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n135 = area1_area3_area5 + area1_area2_area3_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n145 =  area1_area4_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n234 = area2_area3_area4 + area1_area2_area3_area4 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,,
            n235 =  area2_area3_area5 + area1_area2_area3_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n245 = area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n345 = area3_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n1234 = area1_area2_area3_area4 + area1_area2_area3_area4_area5,
            n1235 = area1_area2_area3_area5 + area1_area2_area3_area4_area5,
            n1245 =  area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n1345 = area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n2345 =  area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n12345 =area1_area2_area3_area4_area5,
                category = expnames_tab, 
                euler.d = TRUE, 
                scaled = TRUE, 
                ext.text = FALSE,
                col = col_vec, 
                fill = col_vec,
                cex = rep(2, 31), 
                cat.cex = rep(2.5, 5));
        dev.off();

}


cat("Retrieving list of element per overlap\n");

output_folder_peaks <- paste(output_folder, "peaks_per_circle/", sep="");
checkingOutputFolder(output_folder_peaks);

for(i in 1:length(ol$peaklist)) 
{
	cat("\t", i, "/", length(ol$peaklist), "\n");
	
	if(nchar(names(ol$peaklist)[i]) > 6)
	{
		index_exp_vec <- as.numeric(unlist(lapply(strsplit(unlist(strsplit(names(ol$peaklist)[i],"///")),"peaks"),"[",2)));
		output_file <- paste(expnames_tab[index_exp_vec], collapse="_");
	}
	else
	{
		index_exp <- as.numeric(unlist(lapply(strsplit(names(ol$peaklist)[i],"peaks"),"[",2)));
		output_file <- expnames_tab[index_exp];
	}
	
	gff_table <- data.frame(seqname=as.character(seqnames(ol$peaklist[[i]])), 
			                source="vennDiagram_overlapGFF", 
							feature = as.character(unlist(lapply(elementMetadata(ol$peaklist[[i]])$peakNames, paste,collapse="-"))), 
							start=start(ranges(ol$peaklist[[i]])),
							end=end(ranges(ol$peaklist[[i]])),
							score=0,
							strand=as.character(strand(ol$peaklist[[i]])),
							frame=".",
							group=".")

    write.table(gff_table, file=paste(output_folder_peaks, output_file, ".gff", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
}



cat("Boxplotting occupancy according to different features\n");

TxDb <- makeTxDbFromUCSC(genome= genome_version, tablename="ensGene");
output_folder_occupancy <- paste(output_folder, "genomic_location/", sep="");
checkingOutputFolder(output_folder_occupancy);

#Determining the whole experiment occupancy
for(i in 1:length(gff_GRanges_list)) 
{
	aCR<-assignChromosomeRegion(gff_GRanges_list[[i]], TxDb=TxDb);
	
	if(output_format == "png"){
		png(filename=paste(output_folder_occupancy, expnames_tab[i], "_GeneralgenomicOccupancy.png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder_occupancy, expnames_tab[i], "_GeneralgenomicOccupancy.ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder_occupancy, expnames_tab[i],"_GeneralgenomicOccupancy.pdf"), width=10, height=10)
 }
	barplot(aCR$percentage, main=paste(expnames_tab[i], "_general", sep=""), ylab="percentage");
	dev.off();
}

#Determining occupancy with the categories of the overlap

for(i in 1:length(ol$peaklist)) 
{
	aCR <- assignChromosomeRegion(ol$peaklist[[i]], TxDb=TxDb);
	
	if(nchar(names(ol$peaklist)[i]) > 6)
	{
		index_exp_vec <- as.numeric(unlist(lapply(strsplit(unlist(strsplit(names(ol$peaklist)[i],"///")),"peaks"),"[",2)));
		output_file <- paste(expnames_tab[index_exp_vec], collapse="_");
	}else
	{
		index_exp <- as.numeric(unlist(lapply(strsplit(names(ol$peaklist)[i],"peaks"),"[",2)));
		output_file <- expnames_tab[index_exp];
	}
	
	if(output_format == "png"){
		png(filename=paste(output_folder_occupancy, output_file, "_CirclegenomicOccupancy.png",sep=""), width = 600, height = 600, bg = "transparent")
	}else if(output_format == "ps"){
		cairo_ps(filename=paste(output_folder_occupancy, output_file, "_CirclegenomicOccupancy.ps",sep=""), width = 7, height = 7, bg = "transparent");
	}else{
     pdf(file=paste0(output_folder_occupancy, output_file,"_CirclegenomicOccupancy.pdf"), width=10, height=10)
 }
	barplot(aCR$percentage, main=paste(output_file, "_circle", sep=""), ylab="percentage");
	dev.off();
}


cat("Making heatmaps of the overlapping peaks\n");

if(!is.na(bigwig_vec))
{
	output_folder_heatmaps <- paste(output_folder, "heatmaps/", sep="");
	checkingOutputFolder(output_folder_heatmaps);
	
	for(i in 1:length(ol$peaklist)) 
	{
		cat("\t Generating heatmap for peak combinations ", i, "/", length(ol$peaklist), "\n");
		features <- ol$peaklist[[i]]; 						#Retrieving the last element corresponding to the common peaks between all experiments
		feature.recentered <- feature.center <- features;
		
		if(center_criterion == "coordinates")
		{
			cat("\t\t Defining peak centers according to coordinates\n");
			wid <- width(features); 												#Retrieves the length of all peaks
			start(feature.center) <- start(features) + floor(wid/2); 				#Change start coord to the middle of the peak
		}else #center criterion is max
		{
			cat("\t\t Importing values of reference experiment big wig file\n");
			cvg_refExp <- import(bigwig_vec[grep(ref_exp_max_centering, bigwig_name_vec)], format="BigWig", which=features, as="NumericList");  #lapply(cvglists, function(y){unlist(lapply(y, function(x){return( which(is.na(runValue(x))))}))})
			
			cat("\t\t Defining peak centers according to the maximum value\n");
			index_max <- unlist(lapply(cvg_refExp,which.max));
			start(feature.center) <- start(features) + index_max;	
		}
		
		width(feature.center) <- 1;												#Change the end coord by defining the width to 1
		start(feature.recentered) <- start(feature.center) - upstream_heatmap;
		end(feature.recentered) <- end(feature.center) + downstream_heatmap;
		
		cat("\t\t Importing values of big wig files\n");
		cvglists <- sapply(bigwig_vec, import, format="BigWig", which=feature.recentered, as="RleList");  #lapply(cvglists, function(y){unlist(lapply(y, function(x){return( which(is.na(runValue(x))))}))})
		names(cvglists) <- bigwig_name_vec;
		sig <- featureAlignedSignal(cvglists, feature.center, upstream= upstream_heatmap, downstream= downstream_heatmap);
		
		
		#Replacing negative values by zero: NORMAL IF WARNINGS ARE PRINTED OUT
		sig <- lapply(sig, function(x){x[which(is.na(x))] <- 0; return(x);});
		
		##### Determining name of output file
		
		if(nchar(names(ol$peaklist)[i]) > 6)
		{
			index_exp_vec <- as.numeric(unlist(lapply(strsplit(unlist(strsplit(names(ol$peaklist)[i],"///")),"peaks"),"[",2)));
			output_file <- paste(expnames_tab[index_exp_vec], collapse="_");
		}else
		{
			index_exp <- as.numeric(unlist(lapply(strsplit(names(ol$peaklist)[i],"peaks"),"[",2)));
			output_file <- expnames_tab[index_exp];
		}
		
		cat("\t\t Generating heatmaps\n");
		upper_extreme_vec <- unlist(lapply(sig,max));
		lower_extreme_vec <- unlist(lapply(sig,min));
		
		if(output_format == "png"){
			png(filename=paste(output_folder_occupancy, output_file, "_overlap_heatmap.png",sep=""), width = 600, height = 600, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_occupancy, output_file, "_overlap_heatmap.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste0(output_folder_occupancy, output_file,"_overlap_heatmap.pdf"), width=10, height=10)
  }
		heatmap <- featureAlignedHeatmap(sig, feature.center, upstream= upstream_heatmap, downstream= downstream_heatmap, upper.extreme= upper_extreme_vec, lower.extreme= lower_extreme_vec);
		dev.off();
		
		cat("\t\t Converting heatmap to treeview file\n");
		
		nrow_elements <- unique(as.numeric(unlist(lapply(sig,nrow))));
		
		if(length(nrow_elements) != 1)
		{
			stop("\n Message 1: Problem in the script, this message should not be output, contact the developper\n");
		}
		
		for(i in 1:length(sig))
		{	
			cat("\t\t\t Treating file ", i, length(sig), "\n");
			output_folder_treeview <- paste(output_folder_heatmaps, "treeview/circle_", output_file, "/order_by_", names(sig)[i], "/",  sep="");
			checkingOutputFolder(output_folder_treeview);
			mean_values <- apply(sig[[i]], MARGIN=1, mean);
			indexes_ordered <- order(mean_values);
			
			for(j in 1:length(sig)) 
			{
				ordered_sig <- sig[[j]][indexes_ordered,];
				
				writeTreeviewCluster(fileName = paste(output_folder_treeview, names(sig)[j], "-treeview", sep=""), 
						interpolatedValues = sig[[j]], 
						interpolatedCoordinates = 1:ncol(sig[[j]]), 
						annoNames = 1:nrow(sig[[j]]));
			}
			
		}
		
		cat("\t\t Generating distribution of peaks signal\n");
		if(output_format == "png"){
			png(filename=paste(output_folder_heatmaps, output_file, "_overlap_distribution.png",sep=""), width = 600, height = 600, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(output_folder_heatmaps, output_file, "_overlap_distribution.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste0(output_folder_heatmaps, output_file,"_overlap_distribution.pdf"), width=10, height=10)
  }
		featureAlignedDistribution(sig, feature.center, upstream= upstream_heatmap, downstream= downstream_heatmap, type="l", ylim=c(0, max(unlist(lapply(sig,max))))+1);
		dev.off();
		
		cat("\t\t Retrieving enriched GO terms\n");
		mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= if(organism_name == "mouse") "mmusculus_gene_ensembl" else "hsapiens_gene_ensembl")
		annoData <- getAnnotation(mart);
		annotatedPeak <- annotatePeakInBatch(feature.recentered, AnnotationData = annoData, output="overlapping", PeakLocForDistance="middle", select="all");
		over_condensed <- getEnrichedGO(annotatedPeak, orgAnn= if(organism_name == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db", maxP=0.01, multiAdj=TRUE, minGOterm=1, multiAdjMethod="BH", condense=TRUE);	 
		over_condensed_10 <- getEnrichedGO(annotatedPeak, orgAnn= if(organism_name == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db", maxP=0.01, multiAdj=TRUE, minGOterm=10, multiAdjMethod="BH", condense=TRUE); #Minimum of 10 GO terms to be included in the results
		
		cat("\t\t\t Filtering and sorting by Benjamini-Hochberg p-values\n"); #Normally not necessary but performed as a security
		over_condensed <- BH_filtering_and_sorting(over_condensed);
		over_condensed_10 <- BH_filtering_and_sorting(over_condensed_10);
		
		cat("\t\t Writing genes and GO tables\n");
		
		# Retrieve genes that have overlapping peaks
		ensembl_id <- annotatedPeak$feature[which(!is.na(annotatedPeak$feature))];
		table_gene_to_write <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene", "description"),values= ensembl_id,mart= mart);
		
		output_folder_geneList <- paste(output_folder, "gene_list/", sep="");
		checkingOutputFolder(output_folder_geneList);
		
		write.table(table_gene_to_write, file=paste(output_folder_geneList, output_file, "_geneList.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE);
		
		
		output_folder_GO <- paste(output_folder, "gene_ontology/", sep="");
		checkingOutputFolder(output_folder_GO);
		
		write.table(over_condensed[["bp"]], file=paste(output_folder_GO, output_file, "-biologicalProcess.txt"), sep="\t", quote=FALSE, row.names=FALSE);
		write.table(over_condensed[["mf"]], file=paste(output_folder_GO, output_file, "-molecularFunction.txt"), sep="\t", quote=FALSE, row.names=FALSE);
		write.table(over_condensed[["cc"]], file=paste(output_folder_GO, output_file, "-cellularComponent.txt"), sep="\t", quote=FALSE, row.names=FALSE);
		
		#at least ten peaks per features
		write.table(over_condensed_10[["bp"]], file=paste(output_folder_GO, output_file, "-biologicalProcess10.txt"), sep="\t", quote=FALSE, row.names=FALSE);
		write.table(over_condensed_10[["mf"]], file=paste(output_folder_GO, output_file, "-molecularFunction10.txt"), sep="\t", quote=FALSE, row.names=FALSE);
		write.table(over_condensed_10[["cc"]], file=paste(output_folder_GO, output_file, "-cellularComponent10.txt"), sep="\t", quote=FALSE, row.names=FALSE);
	}
}





