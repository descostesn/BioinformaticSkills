############
# This script performs heatmap of kmean clustering of chip-seq data using annotations.
# Descostes feb 2017
############


library("seqplots");
library("BSgenome");
library("NGSprofiling");
library("Rargs");



################
# PARAMETERS
################


paramsDefinition <- list();


# Required arguments
paramsDefinition[["--bigwigVec"]] <- list(variableName="bigwig_vec", numeric=F, mandatory=T, description="Vector of bigwig files from which the different profiles will be plotted.");
paramsDefinition[["--bigwigNameVec"]] <- list(variableName="bigwig_name_vec", numeric=F, mandatory=T, description="Vector of names corresponding to the bigwig files.");
paramsDefinition[["--bedVec"]] <- list(variableName="bed_vec", numeric=F, mandatory=T, description="Vector of bed files defining the different categories of expression.");
paramsDefinition[["--bedNameVec"]] <- list(variableName="bed_name_vec", numeric=F, mandatory=T, description="Vector of corresponding bed file names.");
paramsDefinition[["--organism"]] <- list(variableName="organism", numeric=F, mandatory=T, description="Should be human, mouse or drosophila.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="Should be hg19, mm10 or dm3.");
paramsDefinition[["--binSize"]] <- list(variableName="bin_size", numeric=T, mandatory=T, description="Binning of the data used on for profiling.");
paramsDefinition[["--profileLengthBefore"]] <- list(variableName="profile_length_before", numeric=T, mandatory=T, description="Positive integer giving the binning uptream the feature.");
paramsDefinition[["--profileLengthAfter"]] <- list(variableName="profile_length_after", numeric=T, mandatory=T, description="Positive integer giving the binning downtream the feature..");
paramsDefinition[["--typeValue"]] <- list(variableName="type_value", numeric=F, mandatory=T, description="pf = start of the feature, mf = midpoint feature, ef= end point feature, af = anchor feature (composite profile).");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Single path to the output folder.");
paramsDefinition[["--nbOfGroups"]] <- list(variableName="nb_of_groups", numeric=T, mandatory=T, description="Single value giving the number of groups for the kmean clustering.");
paramsDefinition[["--includeExpVec"]] <- list(variableName="include_exp_vec", numeric=F, mandatory=T, description="Vector of logical indicating which exp should be used for clustering.", postConversion=as.logical);
paramsDefinition[["--clusteringMethod"]] <- list(variableName="clustering_method", numeric=F, mandatory=T, description="Can be kmeans, hclust or none.");
paramsDefinition[["--sortRows"]] <- list(variableName="sort_rows", numeric=F, mandatory=F, description="Can be FALSE, increasing or decreasing.");
paramsDefinition[["--autoScale"]] <- list(variableName="auto_scale", numeric=F, mandatory=T, description="Logical indicating if colors should be scaled within experiments and not across.", postConversion=as.logical);
paramsDefinition[["--zmin"]] <- list(variableName="zmin", numeric=T, mandatory=T, description="global minimum value on colour key, ignored if ‘autoscale’ is TRUE");
paramsDefinition[["--zmax"]] <- list(variableName="zmax", numeric=T, mandatory=T, description="global maximum value on colour key, ignored if ‘autoscale’ is TRUE");

#optional
paramsDefinition[["--nbPointInterpolation"]] <- list(variableName="nb_point_interpolation", numeric=T, mandatory=F, description="Nb of points used to interpolate the data between the features if type af is used.", default=1000);
paramsDefinition[["--meanOrMedian"]] <- list(variableName="mean_or_median", numeric=F, mandatory=F, description="Mean or median can be used to plot the profile.", default="mean");
paramsDefinition[["--plotScale"]] <- list(variableName="plot_scale", numeric=F, mandatory=F, description="Can be no, linear, log2 or zscore.", default="no");
paramsDefinition[["--colValue"]] <- list(variableName="col_value", numeric=F, mandatory=F, description="Single string giving the color of the gradiant of heatmaps. default=blue", default="blue");





#bigwig_vec <- c("/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0_PolII_mm10_spikedESC.bw",
#                "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-SPT16_mm10_spikedESC.bw",
#                "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5_PolII_mm10_spikedESC.bw", 
#                "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5_SPT16_mm10_spikedESC.bw",
#                "/ifs/home/descon01/data/gary_paper_data/bigwigs/ESC_mouse/spikedin/not_used/WIGfs_Day13_RNAPII_mm10_unireads_elManual176_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#                "/ifs/home/descon01/data/gary_paper_data/bigwigs/ESC_mouse/spikedin/not_used/WIGfs_Day13-SPT16_mm10_unireads_elEst141_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin-save.bw")
#bigwig_name_vec <- c("Day0_RNAPII", "Day0-SPT16", "Day5_RNAPII", "Day5-SPT16", "Day13_RNAPII", "Day13-SPT16");
#bed_vec <- c("/ifs/home/descon01/cluster/Annotations/mouse/mm10/gff_files/allrefseq.gff");
#bed_name_vec  <- c("allrefseq")
#organism <- "mouse"
#genome_version <- "mm10"
#bin_size <- 50
#profile_length_before <- c(10000) 
#profile_length_after <- c(10000)
#type_value <- "pf"
#output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/clustering/decreasing_PolII/ESC/basedallTSS/spiked/PolIISPT16only/");
#nb_of_groups <- c(3)
#include_exp_vec <- c(T, F, F, F, F, F);
#clustering_method <- c("none")
#sort_rows <- "decreasing"
#nb_point_interpolation <- 1000
#mean_or_median <- "mean"
#plot_scale <- "no"
#col_value <- "blue"
#auto_scale <- FALSE


#bigwig_vec <- c("/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0_PolII_mm10_spikedESC.bw",
#        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-SPT16_mm10_spikedESC.bw")
#bigwig_name_vec <- c("Day0_RNAPII", "Day0-SPT16");
#bed_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_macs2/0.04/no_model_broad/Day0_RNAPII_peaks_broadPeak.gff");
#bed_name_vec  <- c("polIIpeaks")
#organism <- "mouse"
#genome_version <- "mm10"
#bin_size <- 50
#profile_length_before <- c(10000) 
#profile_length_after <- c(10000)
#type_value <- "pf"
#output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/clustering/decreasing_PolII/ESC/basedallTSS/spiked/PolIISPT16only/");
#nb_of_groups <- c(3)
#include_exp_vec <- c(T, F);
#clustering_method <- c("none")
#sort_rows <- "decreasing"
#nb_point_interpolation <- 1000
#mean_or_median <- "mean"
#plot_scale <- "no"
#col_value <- "blue"
#auto_scale <- FALSE


################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


if(sort_rows != "increasing" && sort_rows != "decreasing" && sort_rows != "FALSE")
{
	stop("sort_rows should be increasing, decreasing or FALSE\n");
}

if(organism != "human" && organism != "mouse" && organism != "drosophila")
{
	stop("Only human, mouse and drosophila are currently supported\n\n");
}

if(genome_version != "hg19" && genome_version != "mm10" && genome_version != "dm3")
{
	stop("hg19, mm10 and dm3 are currently supported\n");
}

if(type_value != "pf" && type_value != "mf" && type_value != "ef" && type_value != "af")
{
	stop("\n Type value should be pf, mf, ef or af only\n\n");
}


if(length(bigwig_vec) != length(bigwig_name_vec))
{
	stop("\n one name per bigwig should be given\n");
}


if(length(bed_vec) != length(bed_name_vec))
{
	stop("\n one name per bed file should be given\n");
}

if(length(bed_vec) != 1)
{
	stop("Only one bed file can be given as input\n");
}


if(organism == "human")
{
	library(BSgenome.Hsapiens.UCSC.hg19);
}else if(organism == "mouse")
{
	library(BSgenome.Mmusculus.UCSC.mm10)
}else{
	library(BSgenome.Dmelanogaster.UCSC.dm3);
}


if(sort_rows == "FALSE")
{
	sort_rows <- FALSE;
}

checkingOutputFolder(output_folder);

set.seed(45)

#Creating the Setarray object

cat("Creating initial object\n");

complete_array <- getPlotSetArray(tracks = bigwig_vec, 
		features = bed_vec, 
		refgenome = genome_version, 
		bin = bin_size, 
		xmin = profile_length_before, 
		xmax = profile_length_after, 
		xanchored = nb_point_interpolation,
		type = type_value, 
		add_heatmap = TRUE, 
		stat = mean_or_median);


# Plotting all marks by feature

cat("Plotting the heatmap\n");

exp_sorting <- if(all(include_exp_vec)) "all" else paste(bigwig_name_vec[include_exp_vec],collapse="-");
output_file_name <- paste("heatmap_group", nb_of_groups, "_", clustering_method, "_sort", as.character(sort_rows), "_autoscale", as.character(auto_scale), "_sortby", exp_sorting, sep="");

#cairo_ps(filename=paste(output_folder, output_file_name, ".ps",sep=""), width = 15, height = 15, bg = "transparent");
#png(filename=paste(output_folder, output_file_name,".png",sep=""), width = 2000, height = 2000, bg = "transparent")


pdf(file= paste(output_folder, output_file_name, ".pdf", sep=""), onefile=TRUE, paper="a4", colormodel="cmyk", title="seqplots kmean", height=740,width=740)
par(mar=c(1,1,1,1), lwd = 0.2)
coordinates_df <- plotHeatmap(complete_array, 
		labels = bigwig_name_vec,
		legend = TRUE, 
		plotScale = plot_scale, 
		sortrows = sort_rows,
		clusters = nb_of_groups, 
		clstmethod = clustering_method, 
		include = include_exp_vec, 
		autoscale = auto_scale,
  zmin = zmin,
  zmax = zmax,
		ln.v = FALSE, 
		colvec = if(length(col_value) == 1) rep(col_value, length(bigwig_name_vec)) else col_value,
		raster = TRUE,
  ggplot = FALSE,
  clspace = c("white", "blue"),
  cex.lab=2,
cex.axis=1,
cex.legends=2);
dev.off();




# Ouput info about the clustering and make java treeview compatible files

cat("Output gff file of coordinates\n");

gff_to_write <- data.frame(seqname = coordinates_df$chromosome, source="kmean_seqplots", feature = if(!is.null(coordinates_df$metadata_name)) coordinates_df$metadata_name else coordinates_df$metadata_type, start = coordinates_df$start, end = coordinates_df$end, score = coordinates_df$metadata_score, strand = coordinates_df$strand, frame=".", group=coordinates_df$ClusterID)

write.table(gff_to_write, file=paste(output_folder, output_file_name, ".gff", sep=""), sep="\t", quote=F,row.names=F,col.names=F);
