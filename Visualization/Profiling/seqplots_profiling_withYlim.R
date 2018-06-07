############
# This script performs profiling of chip-seq data using annotations.
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
paramsDefinition[["--organism"]] <- list(variableName="organism", numeric=F, mandatory=T, description="Should be human, mouse, drosophila or ant.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="Should be hg19, mm10, dm3 or hsal35.");
paramsDefinition[["--binSize"]] <- list(variableName="bin_size", numeric=T, mandatory=T, description="Binning of the data used on for profiling.");
paramsDefinition[["--profileLengthBefore"]] <- list(variableName="profile_length_before", numeric=T, mandatory=T, description="Positive integer giving the binning uptream the feature.");
paramsDefinition[["--profileLengthAfter"]] <- list(variableName="profile_length_after", numeric=T, mandatory=T, description="Positive integer giving the binning downtream the feature..");
paramsDefinition[["--typeValue"]] <- list(variableName="type_value", numeric=F, mandatory=T, description="pf = start of the feature, mf = midpoint feature, ef= end point feature, af = anchor feature (composite profile).");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Single path to the output folder.");
paramsDefinition[["--scaleBetween01"]] <- list(variableName="scale_between0_1", numeric=F, mandatory=T, description="A boolean defining if values should be scaled between 0 and 1.", postConversion=as.logical);
paramsDefinition[["--colVec_marks"]] <- list(variableName="col_vec_marks", numeric=F, mandatory=T, description="Vector of colors for each experiment.");
paramsDefinition[["--colVec_categories"]] <- list(variableName="col_vec_categories", numeric=F, mandatory=T, description="Vector of colors for each category.");
paramsDefinition[["--outputFormat"]] <- list(variableName="output_format", numeric=F, mandatory=T, description="output format which can be png or ps.");
paramsDefinition[["--nbPointInterpolation"]] <- list(variableName="nb_point_interpolation", numeric=T, mandatory=T, description="Nb of points used to interpolate the data between the features if type af is used. default=10000");
paramsDefinition[["--meanOrMedian"]] <- list(variableName="mean_or_median", numeric=F, mandatory=T, description="Mean or median can be used to plot the profile.default=mean");
paramsDefinition[["--errorEstimates"]] <- list(variableName="error_estimates", numeric=F, mandatory=T, description="A boolean defining if error estimates should be indicated (default: TRUE)", postConversion=as.logical);
#optional
paramsDefinition[["--yLim"]] <- list(variableName="y_lim", numeric=T, mandatory=F, description="Y lim vector of the form c(y1,y2).", default=NULL);


#bigwig_vec <- c("/ifs/home/descon01/data/data_june2017/bam_files/comzit_chipseq/Result_pasha_unireads_SE/AllReads/WIGfs_H3K27me3_G_unireads_elManual151_AThr1_bin50_Scaled_BGSub.bw",
#        "/ifs/home/descon01/data/data_june2017/bam_files/comzit_chipseq/Result_pasha_unireads_SE/AllReads/WIGfs_H3K27me3_W_unireads_elManual151_AThr3_bin50_Scaled_BGSub.bw",
#        "/ifs/home/descon01/data/data_june2017/bam_files/comzit_chipseq/Result_pasha_unireads_SE/AllReads/WIGfs_H3K4me3_G_unireads_elManual121_AThr2_bin50_Scaled_BGSub.bw",
#        "/ifs/home/descon01/data/data_june2017/bam_files/comzit_chipseq/Result_pasha_unireads_SE/AllReads/WIGfs_H3K4me3_W_unireads_elManual141_AThr1_bin50_Scaled_BGSub.bw")
#bigwig_name_vec <- c("H3K27me3_G", "H3K27me3_W", "H3K4me3_G", "H3K4me3_W");
#bed_vec <- c("/ifs/home/descon01/cluster/Annotations/ants/Hsalv3_5/gff_files/hsal_mrna_withchr.gff")
#bed_name_vec <- c("mrna")
#organism <- "ant"
#genome_version <- "HSaltator"
#bin_size <- 50
#profile_length_before <- c(2000)
#profile_length_after <- c(2000)
#type_value <- "af";
#output_folder <- c("/ifs/home/descon01/analysis/comzit_ants/profiles/all_mrna/");
#scale_between0_1 <- "FALSE";
#nb_point_interpolation <- 10000
#mean_or_median <- "mean"
#col_vec_marks <- c("#7570B3", "#E7298A", "#66A61E", "#E6AB02");
#col_vec_categories <- c("#7570B3"); 
#output_format <- "png";


#bigwig_vec <-c("/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_Flag-hdgf2_293_native_filtered_unireads_elEst121_AThr7_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_Flag-ledgf_293_native_filtered_unireads_elEst131_AThr6_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K27me3_filtered_unireads_elEst191_AThr4_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K36me2_filtered_unireads_elManual184_AThr4_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K36me3_filtered_unireads_elManual173_AThr5_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_merge_spt16-1-2_filtered_filtered_unireads_elManual162_AThr4_bin50_Scaled_BGSub.bw",
#				"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_RNAPolII_filtered_unireads_elEst171_AThr3_bin50_Scaled_BGSub.bw")
#
#bigwig_name_vec <- c("Hdgf", "Ledgf", "H3K27me3", "H3K36me2", "H3K36me3", "spt16", "PolII");
#
#bed_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category1/refseq.gff", 
#		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category2/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category3/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category4/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category5/refseq.gff",
#		"/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/overalp_categories_refseq/heatmapIntervals_vs_refseq/two_experiments/category6/refseq.gff")
#
#bed_name_vec <- c("category1-genes", "category2-genes", "category3-genes", "category4-genes", "category5-genes", "category6-genes")
#
#organism <- "human"
#genome_version <- "hg19"
#bin_size <- 50
#profile_length_before <- c(5000)
#profile_length_after <- c(5000)
#type_value <- "af";
#output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/clustering/kmean/293/no_model_broad_PolII_peaks_0.03/gary/10000/profiles/per_category/all_refseq_in_intervals/gary/5000/scale_between0_1/");
#scale_between0_1 <- TRUE;
#nb_point_interpolation <- 10000
#mean_or_median <- "mean"

#bigwig_vec <- c("/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_Flag-hdgf2_293_native_filtered_unireads_elEst121_AThr7_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_Flag-ledgf_293_native_filtered_unireads_elEst131_AThr6_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K27me3_filtered_unireads_elEst191_AThr4_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K36me2_filtered_unireads_elManual184_AThr4_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_H3K36me3_filtered_unireads_elManual173_AThr5_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_merge_spt16-1-2_filtered_filtered_unireads_elManual162_AThr4_bin50_Scaled_BGSub.bw",
#		"/ifs/home/descon01/data/gary_paper_data/bigwigs/293TRex/WIGfs_RNAPolII_filtered_unireads_elEst171_AThr3_bin50_Scaled_BGSub.bw");
#bigwig_name_vec <- c("Hdgf", "Ledgf", "H3K27me3", "H3K36me2", "H3K36me3", "spt16", "PolII");
#bed_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/rnaseq/order_by_expression/different_classification/rep1mergetech/3-groups/edgeR_rpkm/fisher/fisher-category1.bed",
#		"/ifs/home/descon01/analysis/fact_ledgf/rnaseq/order_by_expression/different_classification/rep1mergetech/3-groups/edgeR_rpkm/fisher/fisher-category2.bed",
#		"/ifs/home/descon01/analysis/fact_ledgf/rnaseq/order_by_expression/different_classification/rep1mergetech/3-groups/edgeR_rpkm/fisher/fisher-category3.bed");
#bed_name_vec <- c("category-1", "category-2", "category-3")
#organism <- "human";
#genome_version <- "hg19";
#bin_size <- 50;
#profile_length_before <- 1000;
#profile_length_after <- 1000;
#type_value <- "af"; #pf = start of the feature, mf = midpoint feature, ef= end point feature, af = anchor feature (composite profile) 
#output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/with_code/rep1mergetech/3-groups/edgeR_rpkm/fisher/"
###only used if type is "af"
#nb_point_interpolation <- 10000;
#mean_or_median <- "mean";

################


################
# FUNCTION
################


range01 <- function(x){(x-min(x))/(max(x)-min(x))};


################


##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


if(organism != "human" && organism != "mouse" && organism != "drosophila" && organism != "ant")
{
    stop("Only human, mouse, drosophila and ant are currently supported\n\n");
}

if(genome_version != "hg19" && genome_version != "mm10" && genome_version != "dm3" && genome_version != "hsal35")
{
    stop("hg19, mm10, dm3 and hsal35 are currently supported\n");
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


if(organism == "human")
{
    library(BSgenome.Hsapiens.UCSC.hg19);
}else if(organism == "mouse")
{
    library(BSgenome.Mmusculus.UCSC.mm10)
}else if (organism == "drosophila"){
    library(BSgenome.Dmelanogaster.UCSC.dm3);
}else{
    library(BSgenome.HSaltator.Lab.hsal35);
}


if(!is.null(col_vec_marks) && length(col_vec_marks) != length(bigwig_vec))
{
    stop("One color should be attributed to each experiment\n");
}

if(!is.null(col_vec_categories) && length(col_vec_categories) != length(bed_vec))
{
    stop("One color should be attributed to each category\n");
}


if(output_format != "png" && output_format != "ps" && output_format != "pdf")
{
    stop("output_format should be png, ps or pdf\n");
}

checkingOutputFolder(output_folder);

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
        add_heatmap = FALSE, 
        stat = mean_or_median);




if(scale_between0_1)
{
    cat("Scaling values between 0 and 1\n");
    for(i in 1:length(complete_array$data))
    {
        for(j in  1:length(complete_array$data[[i]]))
        {
            values_mean <- complete_array$data[[i]][[j]]$means;
            co <- diff(range(values_mean));
            complete_array$data[[i]][[j]]$means <- range01(values_mean);
            co <- diff(range(complete_array$data[[i]][[j]]$means))/co; 
            complete_array$data[[i]][[j]]$stderror <- co*complete_array$data[[i]][[j]]$stderror;
            complete_array$data[[i]][[j]]$conint <- co*complete_array$data[[i]][[j]]$conint;
        }
    }
}

# Plotting all marks by feature

cat("Plotting all marks by feature\n");


for(i in 1:length(bed_vec)) 
{
    cat("\t Plotting feature ", i, "/", length(bed_vec), "\n");
    current_selection <- complete_array$get(i, 1:length(bigwig_vec));
    
    if(output_format == "png"){
        png(filename=paste(output_folder, bed_name_vec[i], "-", type_value, "-", paste(bigwig_name_vec, collapse="_"), "-allMarks.png",sep=""), width = 600, height = 600, bg = "transparent");
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, bed_name_vec[i], type_value, "-", paste(bigwig_name_vec, collapse="_"), "-allMarks.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste(output_folder, bed_name_vec[i], type_value, "-", paste(bigwig_name_vec, collapse="_"), "-allMarks.pdf",sep=""), width=10, height=10)
    } 
    plotAverage(current_selection, 
            keepratio = FALSE,
            labels = bigwig_name_vec,  
            plotScale = "linear", 
            type = "full",
            cex.lab = 9,
            colvec = col_vec_marks,
            error.estimates = error_estimates,
            ylim = y_lim);
    dev.off();
}



# Plotting all features by marks

cat("Plotting all features by marks\n");

for(i in 1:length(bigwig_name_vec)) 
{
    cat("\t Plotting feature ", i, "/", length(bigwig_name_vec), "\n");
    current_selection <- complete_array$get(1:length(bed_vec),i);
    
    if(output_format == "png"){
        png(filename=paste(output_folder, bigwig_name_vec[i], "-", type_value, "-allCategories.png",sep=""), width = 600, height = 600, bg = "transparent");
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, bigwig_name_vec[i], "-", type_value, "-allCategories.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste(output_folder, bigwig_name_vec[i], "-", type_value, "-allCategories.pdf",sep=""), width=10, height=10)
    } 
    plotAverage(current_selection, 
            keepratio = FALSE,
            labels = bed_name_vec,  
            plotScale = "linear", 
            type = "full",
            cex.lab = 9,
            colvec = col_vec_categories,
            error.estimates = error_estimates,
            ylim = y_lim);
    dev.off();
}


#Plotting each mark with each feature

cat("Plotting each mark with each feature\n");

for(i in 1:length(bed_name_vec)) 
{
    for(j in 1:length(bigwig_name_vec)) 
    {
        current_selection <- complete_array$get(i, j);
        
        cat("\t Plotting: ", bigwig_name_vec[j], "/", bed_name_vec[i], "\n");
        
        if(output_format == "png"){
            png(filename=paste(output_folder, bigwig_name_vec[j], "-", bed_name_vec[i], "-", type_value, ".png",sep=""), width = 600, height = 600, bg = "transparent");
        }else if(output_format == "ps"){
            cairo_ps(filename=paste(output_folder, bigwig_name_vec[j], "-", bed_name_vec[i], "-", type_value, ".ps",sep=""), width = 7, height = 7, bg = "transparent");
        }else{
            pdf(file=paste(output_folder, bigwig_name_vec[j], "-", bed_name_vec[i], "-", type_value, ".pdf",sep=""), width=10, height=10)
        } 
        
        plotAverage(current_selection, 
                keepratio = FALSE,
                labels = bigwig_name_vec[j],  
                plotScale = "linear", 
                type = "full",
                cex.lab = 9,
                colvec = col_vec_marks[j],
                error.estimates = error_estimates,
                ylim = y_lim);
        dev.off();
    }
}
