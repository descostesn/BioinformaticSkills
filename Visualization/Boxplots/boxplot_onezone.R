###########
# This script takes a GFF/BED and bigwigs as input and generate a boxplot of mean values on the given intervals.
# Descostes Feb 2018
###########

library(rtracklayer);
library(GenomicRanges);
library("Rargs");


################
# PARAMETERS
################



#parameters defined from the command line using RIO
paramsDefinition <- list();

# Required arguments
paramsDefinition[["--bigwigsVec"]] <- list(variableName="bigwigs_vec", numeric=F, mandatory=T, description="A vector of bigwig file pathes.");
paramsDefinition[["--expnamesVec"]] <- list(variableName="expname_vec", numeric=F, mandatory=T, description="A vector of experiments name corresponding to the bigwig files.");
paramsDefinition[["--expnamesRatioVec"]] <- list(variableName="expname_ratio_vec", numeric=F, mandatory=T, description="A vector of labels for the log2(FC) boxplots.");
paramsDefinition[["--colorsTab"]] <- list(variableName="colors_tab", numeric=F, mandatory=T, description="Vector of colors for each bigwig file.");
paramsDefinition[["--mainTitle"]] <- list(variableName="main_title", numeric=F, mandatory=T, description="Title printed on top of boxplot.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Unique string giving the path to the output folder.");
paramsDefinition[["--yLabelVec"]] <- list(variableName="y_label_vec", numeric=F, mandatory=T, description="Title of the y axis.");
paramsDefinition[["--yLabelRatioVec"]] <- list(variableName="y_label_ratio_vec", numeric=F, mandatory=T, description="Title of the y axis for the log2(FC) plot.");
paramsDefinition[["--gffVec"]] <- list(variableName="gff_vec", numeric=F, mandatory=T, description="Intervals to extract the bigwig values from.");
paramsDefinition[["--nameOutputFileNoextVec"]] <- list(variableName="name_outputFile_noext_vec", numeric=F, mandatory=T, description="Output file name without extension.");



#bigwigs_vec <- c("/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/0h_K27M_WT/WIGfs_K36me2_K27WT0Hr_filtered_unireads_elManual211_AThr2_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_WT/WIGfs_H3K36me2_6hr-H3WT-293TREX_filtered_unireads_elManual220_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_WT/WIGfs_H3K36me2_12hr-H3WT-293TREX_filtered_unireads_elManual220_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_WT/WIGfs_H3K36me2_24hr-H3WT-293TREX_filtered_unireads_elManual220_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/0h_K27M_WT/WIGfs_K36me2_K27M0Hr_filtered_unireads_elManual211_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_M/WIGfs_H3K36me2_6hr-H3K27M-293TREX_filtered_unireads_elManual211_AThr2_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_M/WIGfs_H3K36me2_12hr-H3K27M-293TREX_filtered_unireads_elManual211_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
#        "/home/descostes/Documents/analysis/K27M_project/bigwigs/feb2018/timecourse_M/WIGfs_H3K36me2_24hr-H3K27M-293TREX_filtered_unireads_elManual211_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw")


#bigwigs_vec <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/WT_EED.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/F97A_EED.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/Y358A_EED.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/Y365A_EED.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/KO_EED.bw")

#bigwigs_vec <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/WT_H3K27me3.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/F97A_H3K27me3.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/Y358A_H3K27me3.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/Y365A_H3K27me3.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/KO_H3K27me3.bw")
#bigwigs_vec <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/SUZ12_E14_WT.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/SUZ12_E14_Y365A.bw",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/bw_files/SUZ12_E14_KO.bw")


#expname_vec <- c("WT_EED",
#        "F97A_EED",
#        "Y358A_EED",
#        "Y365A_EED",
#        "KO_EED")
#expname_vec <- c("WT_H3K27me3",
#        "F97A_H3K27me3",
#        "Y358A_H3K27me3",
#        "Y365A_H3K27me3",
#        "KO_H3K27me3")
#expname_vec <- c("SUZ12_E14_WT",
#        "SUZ12_E14_Y365A",
#        "SUZ12_E14_KO")


#plot_ratio <- FALSE
#expname_ratio_vec <- NA
#
#colors_tab <- c(rgb( 0,0,139, maxColorValue = 255),
#        rgb( 0,0,0, maxColorValue = 255),
#        rgb( 128,0,0, maxColorValue = 255),
#        rgb( 0,100,0, maxColorValue = 255),
#        rgb( 75,0,130, maxColorValue = 255))
#colors_tab <- c(rgb( 0,0,139, maxColorValue = 255),
#        rgb( 0,0,0, maxColorValue = 255),
#        rgb( 128,0,0, maxColorValue = 255),
#        rgb( 0,100,0, maxColorValue = 255),
#        rgb( 75,0,130, maxColorValue = 255))
#colors_tab <- c(rgb( 0,0,139, maxColorValue = 255),
#        rgb( 0,100,0, maxColorValue = 255),
#        rgb( 75,0,130, maxColorValue = 255))



#output_folder <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/results/EED/")
#output_folder <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/results/K27me3/")
#output_folder <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/results/SUZ12/")


#y_label_vec <- c("Average Binding on ")
#
#gff_vec <- c("/home/descostes/Documents/analysis/EED-mutant/boxplots/annotations/ALL_peaks.gff",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/annotations/GroupIII_spreadingPeaks.gff",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/annotations/GroupII_weakpeaks.gff",
#        "/home/descostes/Documents/analysis/EED-mutant/boxplots/annotations/GroupI_strongpeaks.gff")
#
#main_title_vec <- c("ALL peaks",
#        "Spreading Peaks",
#        "Weak nucleation sites",
#        "Strong nucleation sites")
#
#name_outputFile_noext_vec <- c("ALL_peaks",
#        "spreadingPeaks",
#        "weaknucleationsites",
#        "strongnucleationsites")

#y_label_ratio_vec <- c("Log2(FC)")



################

################
# FUNCTION
################


retrieveMatBindingList <- function(files_vec, gff_vec, interpolation_number, 
        expname_vec, ignore_strand, verbose){
    
    current_gff <- read.table(gff_vec, stringsAsFactors=FALSE)
    current_gff_GRanges <- GRanges(seqnames=current_gff[,1], 
            ranges=IRanges(start= current_gff[,4], 
                    end = current_gff[,5], 
                    names = current_gff[,2]), 
            strand=current_gff[,7])
    
    mat_list <- mapply(function(file_path, expname){
                
                if(verbose)
                    message("\t\t Processing ", expname)
                
                mat_tmp <- summary(BigWigFile(file_path), 
                        which = current_gff_GRanges, 
                        size = interpolation_number, 
                        type = "mean", as = "matrix")
                
                if (!ignore_strand)
                    mat_tmp[
                            as.character(strand(current_gff_GRanges))=='-',] <-
                            mat_tmp[
                                    as.character(strand(current_gff_GRanges))==
                                            '-', ncol(mat_tmp):1]
                return(mat_tmp)
            }, files_vec, expname_vec, SIMPLIFY = FALSE)
    
    names(mat_list) <- expname_vec
    return(mat_list)
}



################


##############
# MAIN
##############


# Retreives the parameters
#getParams(paramsDefinition);


if(length(bigwigs_vec) != length(expname_vec))
    stop("one exp name should be given per exp")

if(length(gff_vec) != length(name_outputFile_noext_vec))
    stop("one output file name is needed per gff file")

if(!file.exists(output_folder))
{
	dir.create(output_folder, recursive = TRUE)
}


## Retrieving the values from bigwigs on gff intervals

for(i in 1:length(gff_vec)) 
{
    cat("Annotation ", i, "/", length(gff_vec), "\n")
    current_gff <- gff_vec[i];
    name_outputFile <- name_outputFile_noext_vec[i]
    current_gff_loaded <- read.table(current_gff, stringsAsFactors=FALSE)
    main_title <- main_title_vec[i]
    
    zz <- file(paste(output_folder, name_outputFile, "-", main_title, ".log",sep=""), open="wt");
    sink(zz);
    sink(type="message");
    
    cat("\t Retrieving matrices\n");
    mat_list <- retrieveMatBindingList(bigwigs_vec, current_gff, 10000, expname_vec, ignore_strand = TRUE, verbose = TRUE)
    
    cat("\t Computing mean values\n");
    mean_list <- lapply(mat_list, function(mat){return(apply(mat, MARGIN=1, mean))})
    rm(mat_list)
    gc()
    
    cat("\t Building matrix of mean values\n")
    result_mat <- do.call(cbind, mean_list)
    nb_row_before <- nrow(result_mat)
    result_mat <- na.omit(result_mat)
    nb_rows_removed_na <- abs(nrow(result_mat)-nb_row_before)
    cat("\n\t\t############# Nb rows removed due to NA: ", nb_rows_removed_na, "\n")
    
    
    
    cat("\t Plotting boxplots separately\n")
    if(is.null(colors_tab))
        colors_tab <- c("darkblue", "darkgreen", "darkred", "darkmagenta",
                "darkgray","darkorange", "darkcyan", "black",
                rainbow(length(bigwigs_vec)-8))
    
    png(filename=paste(output_folder, name_outputFile,".png",sep=""), width = 600, height = 600, bg = "transparent")
    boxplot(result_mat, main = main_title, ylab= paste(y_label_vec, strsplit(basename(current_gff), "\\.")[[1]][1], sep=""), col= colors_tab, outline=TRUE, names = expname_vec, las=2, cex.axis=0.7, notch=TRUE);
    dev.off()
    
    png(filename=paste(output_folder, name_outputFile,"_nooutliers.png",sep=""), width = 600, height = 600, bg = "transparent")
    boxplot(result_mat, main = main_title, ylab= paste(y_label_vec, strsplit(basename(current_gff), "\\.")[[1]][1], sep=""), col= colors_tab, outline=FALSE, names = expname_vec, las=2, cex.axis=0.7, notch=TRUE);
    dev.off()
    
    write.table(result_mat, file=paste0(output_folder, name_outputFile, "_matrix.txt"), sep="\t", quote=F, col.names=F,row.names=F)
    
    if(plot_ratio){
        cat("\t Plotting ratio\n")
        index_1 <- seq(from=1, to=length(mean_list), by=2)
        index_2 <- seq(from=2, to=length(mean_list), by=2)
        
        if(length(index_1) != length(index_2))
            stop("pb building index1 and 2.")
        
        cat("\t\t Filtering zero from matrix\n")
        denominator_mat <- result_mat[,index_1]
        keep_row_zero <- apply(denominator_mat, MARGIN=1, function(x){ all(x != 0)})
        nb_rows_removed_zero <- length(which(!keep_row_zero))
        cat("\t\t\t########## Number of rows removed because of zero: ", nb_rows_removed_zero, "\n")
        result_mat <- result_mat[keep_row_zero,]
        
        result_mat_ratio_list <- list();
        
        for(j in 1:length(index_1)) 
        {
            mean_values1 <- result_mat[,index_1[j]]
            mean_values2 <- result_mat[,index_2[j]]
            
            result_mat_ratio_list[[j]] <- log2(mean_values2/mean_values1)
            
            nb_rows_removed_na_log2 <- length(which(is.na(result_mat_ratio_list[[j]])))
            cat("\t\t ########## Number of points discarded due to a negative ratio for boxplot", j, ":", nb_rows_removed_na_log2, "/", length(result_mat_ratio_list[[j]]), "\n")
        }
        
        result_mat_ratio <- do.call(cbind, result_mat_ratio_list)
        colors_tab_ratio <- c("darkblue", "darkgreen", "darkred", "darkmagenta",
                "darkgray","darkorange", "darkcyan", "black", rainbow(length(index_1)-8))
        
        cat("\t\t Generating boxplots\n")
        png(filename=paste(output_folder, name_outputFile,"ratio-nooutliers.png",sep=""), width = 600, height = 600, bg = "transparent")
        boxplot(result_mat_ratio, main = main_title, ylab= y_label_ratio_vec, col = colors_tab_ratio, outline=FALSE, names = expname_ratio_vec, las=2, cex.axis=0.7, notch=TRUE)
        dev.off()
        
        png(filename=paste(output_folder, name_outputFile,"ratio.png",sep=""), width = 600, height = 600, bg = "transparent")
        boxplot(result_mat_ratio, main = main_title, ylab= y_label_ratio_vec, col = colors_tab_ratio, outline=TRUE, names = expname_ratio_vec, las=2, cex.axis=0.7, notch=TRUE)
        dev.off()
        
        
        cat("\n\n\t\t ######## Total number of values discarded during the analysis: ", nb_rows_removed_na + nb_rows_removed_na_log2 + nb_rows_removed_zero, "/", nrow(current_gff_loaded), "\n")
    }else{
        cat("\n\n\t\t ######## Total number of values discarded during the analysis: ", nb_rows_removed_na, "/", nrow(current_gff_loaded), "\n")
    }
    sink()
    rm(mean_list)
    rm(result_mat)
    gc()
}

