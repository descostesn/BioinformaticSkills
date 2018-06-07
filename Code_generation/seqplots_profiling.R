
bigwig_vec <- c(paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K27me3_myoblast_filtered_unireads_elEst191_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K27me3_myotube_filtered_unireads_elManual191_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K36me2_myoblast_filtered_unireads_elManual163_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K36me2_myotube_filtered_unireads_elManual163_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K36me3_myoblast_filtered_unireads_elManual115_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_H3K36me3_myotube_filtered_unireads_elManual115_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_HDGF2_myoblast_filtered_unireads_elEst181_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_HDGF2_myotube_filtered_unireads_elEst181_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_LEDGF_myoblast_filtered_unireads_elManual171_AThr7_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_LEDGF_myotube_filtered_unireads_elManual171_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_RNAPII_myoblast_filtered_unireads_elEst181_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_RNAPII_myotube_filtered_unireads_elEst171_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "),
        paste("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_SPT16_myoblast_filtered_unireads_elEst181_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
                "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/Result_pasha_unireads_SE/AllReads/WIGfs_SPT16_myotube_filtered_unireads_elEst141_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw", sep=" "))


bigwig_name_vec <- c(paste("H3K27me3_myoblast",
                "H3K27me3_myotube", sep=" "),
        paste("H3K36me2_myoblast",
                "H3K36me2_myotube", sep=" "),
        paste("H3K36me3_myoblast",
                "H3K36me3_myotube", sep=" "),
        paste("HDGF2_myoblast",
                "HDGF2_myotube", sep=" "),
        paste("LEDGF_myoblast",
                "LEDGF_myotube", sep=" "),
        paste("RNAPII_myoblast",
                "RNAPII_myotube", sep=" "),
        paste("SPT16_myoblast",
                "SPT16_myotube", sep=" "))

bed_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low20percent/peaks_per_circle/polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low20percent/peaks_per_circle/polII-MB_polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low20percent/peaks_per_circle/polII-MB_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low5percent/peaks_per_circle/polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low5percent/peaks_per_circle/polII-MB_polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/low5percent/peaks_per_circle/polII-MB_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20-80percent/peaks_per_circle/polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20-80percent/peaks_per_circle/polII-MB_polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20-80percent/peaks_per_circle/polII-MB_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20percent/peaks_per_circle/polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20percent/peaks_per_circle/polII-MB_polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top20percent/peaks_per_circle/polII-MB_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top5percent/peaks_per_circle/polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top5percent/peaks_per_circle/polII-MB_polII-MT_original.gff",
        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/fromProfiles_on_genes/myoblast_myotube/based_rnapII_myoblast_or_myotube/polII_comparison/quantile/AVG_LEVEL_TSS/top5percent/peaks_per_circle/polII-MB_original.gff")

bed_name_vec <- c("low20percent_polII-MT",
        "low20percent_polII-MB_polII-MT",
        "low20percent_polII-MB",
        "low5percent_polII-MT",
        "low5percent_polII-MB_polII-MT",
        "low5percent_polII-MB",
        "top20-rcent_polII-MT",
        "top20-80percent_polII-MB_polII-MT",
        "top20-80percent_polII-MB",
        "top20percent_polII-MT",
        "top20percent_polII-MB_polII-MT",
        "top20percent_polII-MB",
        "top5percent_polII-MT",
        "top5percent_polII-MB_polII-MT",
        "top5percent_polII-MB")


organism <- "human"
genome_version <- "hg19"
bin_size <- 50
profile_length_before <- c(2000)
profile_length_after <- c(2000)
type_value <- "af";
output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low20percent/polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low20percent/polII-MB_polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low20percent/polII-MB/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low5percent/polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low5percent/polII-MB_polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/low5percent/polII-MB/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20-80percent/polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20-80percent/polII-MB_polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20-80percent/polII-MB/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20percent/polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20percent/polII-MB_polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top20percent/polII-MB/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top5percent/polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top5percent/polII-MB_polII-MT/",
        "/ifs/home/descon01/analysis/fact_ledgf/profiles/seqplot/myoblast_myotube/based_overlap_PolII_quantile_selection/AVG_LEVEL_TSS/top5percent/polII-MB/")

scale_01 <- "FALSE";
#NULL if no colors

col_vec_marks <- c(paste("blue4", "darkred", sep=" ")); 

col_vec_categories <- c(paste("#D53E4F", sep=" ")); 

output_format <- "pdf";

nb_point_interpolation <- 10000
mean_or_median <- "mean"
error_estimates <- TRUE


to_write <- vector()

for (i in seq_len(length(profile_length_after))) {
    for (j in seq_len(length(bed_vec))) {
        to_write <- c(to_write, paste(bigwig_vec, bigwig_name_vec, bed_vec[j], bed_name_vec[j], organism, genome_version, bin_size, profile_length_before[i], profile_length_after[i], type_value, 
                        paste(output_folder[j], profile_length_after[i], "/", sep=""), scale_01, col_vec_marks, col_vec_categories, output_format, nb_point_interpolation, mean_or_median, error_estimates,
                        sep=";"))
    }
}

message("The number of lines should be ", length(bigwig_vec)*length(profile_length_after)*length(bed_vec))

write(to_write, file="gary_overlapPolII_quantileSelection_allmarkscorresponding-tss.conf", ncolumns=1)

