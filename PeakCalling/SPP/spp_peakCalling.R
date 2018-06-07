#############
# This script performs peak calling with spp.
# Descostes March 2017
#############


library("spp");
library("snow");
library("NGSprofiling");
library("Rargs");

# with r/3.3.2


################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();

paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Number of cpu to use for normalization.");
paramsDefinition[["--bamFile"]] <- list(variableName="bam_file", numeric=F, mandatory=T, description="Single string containing file path to experiment bam file.");
paramsDefinition[["--inputFile"]] <- list(variableName="input_file", numeric=F, mandatory=T, description="Single string containing file path to input bam file.");
paramsDefinition[["--outputFolderForPeaks"]] <- list(variableName="output_folder_for_peaks", numeric=F, mandatory=T, description="Single string giving path to the peaks output folder.");
paramsDefinition[["--outputFolderForWigs"]] <- list(variableName="output_folder_for_wigs", numeric=F, mandatory=T, description="Single string giving path to the wigs output folder.");
paramsDefinition[["--expname"]] <- list(variableName="expname", numeric=F, mandatory=T, description="Single string giving the name of the experiment.");
paramsDefinition[["--minDistanceBetweenPeaks"]] <- list(variableName="min_distance_between_peaks", numeric=T, mandatory=T, description="Minimum distance that should separate two peaks.");

#optional param
paramsDefinition[["--windowSizeForPeakDetection"]] <- list(variableName="window_size_forPeakDetection", numeric=T, mandatory=F, description="Size of the window used for peak detection.", default=1000);
paramsDefinition[["--zscoreThreshold"]] <- list(variableName="zscore_threshold", numeric=T, mandatory=F, description="Z-score threshold for detection.", default = 3);
paramsDefinition[["--fdrValue"]] <- list(variableName="fdr_value", numeric=T, mandatory=F, description="False discovery rate score.", default = 0.01);



#nb_cpu <- 2;
#bam_file <- "/ifs/home/descon01/data/data_october2016/bam_files/james_chipseq/EZH2_H33WT_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam";
#input_file <- "/ifs/home/descon01/data/data_october2016/bam_files/james_chipseq/input_H33WT_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam";
#output_folder_for_peaks <- "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_spp/";
#output_folder_for_wigs <- "/ifs/home/descon01/data/data_october2016/bam_files/james_chipseq/Result_SPP/";
#
#expname <- "EZH2_H33WT";
#min_distance_between_peaks <- 200;
#window_size_forPeakDetection <- 1e3
#zscore_threshold <- 3
#fdr_value <- 1e-2

################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);



checkingOutputFolder(output_folder_for_peaks);
checkingOutputFolder(output_folder_for_wigs);


if(length(bam_file) > 1 || length(input_file) > 1)
{
	stop("\n This script takes a single bam and input file\n");
}


cluster <- makeCluster(nb_cpu);


cat("Reading experiment and input files\n");

chip_data <- read.bam.tags(bam_file);
input_data <- read.bam.tags(input_file);

cat("Computing the binding characteristics\n");

binding_characteristics <- get.binding.characteristics(chip_data, srange=c(36,200),bin=5,cluster=cluster);

# print out binding peak separation distance
print(paste("binding peak separation distance =", binding_characteristics$peak$x));

cat("Plot cross-correlation profile\n");

pdf(file=paste(output_folder_for_wigs, expname, "-crosscorrelation.pdf", sep=""),width=5,height=5);
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot(binding_characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
abline(v=binding_characteristics$peak$x,lty=2,col=2)
dev.off();


cat("Select informative tags based on the binding characteristics\n");

chip_data <- select.informative.tags(chip_data, binding_characteristics);
input_data <- select.informative.tags(input_data, binding_characteristics);


cat("Restrict or remove singular positions with very high tag counts\n");

chip_data <- remove.local.tag.anomalies(chip_data);
input_data <- remove.local.tag.anomalies(input_data);


cat("Outputing density wig file\n");

# output smoothed tag density (subtracting re-scaled input) into a WIG file
# note that the tags are shifted by half of the peak separation distance
tag_shift <- round(binding_characteristics$peak$x/2);
smoothed_density <- get.smoothed.tag.density(chip_data,control.tags=input_data,bandwidth=200,step=100,tag.shift=tag_shift);
writewig(smoothed_density, paste(output_folder_for_wigs, expname, "_density.wig", sep=""), expname);
rm(smoothed_density);
gc(verbose = FALSE);


cat("Outputing fold-enrichment wig file\n");

smoothed_enrichment_estimate <- get.smoothed.enrichment.mle(chip_data, input_data, bandwidth=200, step=100, tag.shift= tag_shift);
writewig(smoothed_enrichment_estimate, paste(output_folder_for_wigs, expname, "_enrichmentMaxLikelihoodLog2.wig", sep=""), expname);


cat("Outputing enrichment estimate wig file\n");

# output conservative enrichment estimates
# alpha specifies significance level at which confidence intervals will be estimated
enrichment_estimates <- get.conservative.fold.enrichment.profile(chip_data, input_data, fws=500,step=100,alpha=0.01);
writewig(enrichment_estimates, paste(output_folder_for_wigs, expname, "_enrichmentEstimates_log2scale.wig", sep=""), expname);
rm(enrichment_estimates);
gc(verbose=FALSE);


cat("Writing broad peaks file\n");

broad_clusters <- get.broad.enrichment.clusters(chip_data,input_data, window.size= window_size_forPeakDetection, z.thr= zscore_threshold,tag.shift=round(binding_characteristics$peak$x/2))
write.broadpeak.info(broad_clusters, paste(output_folder_for_peaks, expname, ".broadPeak", sep=""));


cat("Writing narrow peaks positions\n");

# binding detection parameters
# desired FDR (1%). Alternatively, an E-value can be supplied to the method calls below instead of the fdr parameter
fdr <- fdr_value; 
# the binding.characteristics contains the optimized half-size for binding detection window
detection_window_halfsize <- binding_characteristics$whs;

# determine binding positions using wtd method
bp <- find.binding.positions(signal.data=chip_data,control.data=input_data,fdr=fdr,whs=detection_window_halfsize,cluster=cluster, min.dist = min_distance_between_peaks, enrichment.z = zscore_threshold)

print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));

# output detected binding positions
output.binding.results(bp, paste(output_folder_for_peaks, expname, "_binding_positions.txt", sep=""));

#Alternatively, the binding positions can be determined using MTC method (referred here as lwcc):
bp_MTC <- find.binding.positions(signal.data=chip_data,control.data=input_data,fdr=fdr,method=tag.lwcc,whs=detection_window_halfsize,cluster=cluster, enrichment.z = zscore_threshold)

print(paste("detected",sum(unlist(lapply(bp_MTC$npl,function(d) length(d$x)))),"peaks"));

# output detected binding positions
output.binding.results(bp_MTC, paste(output_folder_for_peaks, expname, "_binding_positions_MTC.txt", sep=""));

bp <- add.broad.peak.regions(chip_data,input_data,bp,window.size= window_size_forPeakDetection, z.thr= zscore_threshold);
bp_MTC <- add.broad.peak.regions(chip_data,input_data,bp_MTC,window.size= window_size_forPeakDetection, z.thr= zscore_threshold);

# output using narrowPeak format
write.narrowpeak.binding(bp, paste(output_folder_for_peaks, expname, ".narrowPeak", sep=""));
write.narrowpeak.binding(bp_MTC, paste(output_folder_for_peaks, expname, "_MTC.narrowPeak", sep=""));














