###############
# This script generates a file containing lines to paste in UCSC genome browser to visualize bigwig.
# Descostes, January 2016, update may 2018
###############



#########################
# PARAMETERS
#########################


bigwig_file_name_vec <- c("WIGfs_H3K27me3_myoblast_filtered_unireads_elEst191_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_H3K27me3_myotube_filtered_unireads_elManual191_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_H3K36me2_myoblast_filtered_unireads_elManual163_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_H3K36me2_myotube_filtered_unireads_elManual163_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_H3K36me3_myoblast_filtered_unireads_elManual115_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_H3K36me3_myotube_filtered_unireads_elManual115_AThr6_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_HDGF2_myoblast_filtered_unireads_elEst181_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_HDGF2_myotube_filtered_unireads_elEst181_AThr4_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_LEDGF_myoblast_filtered_unireads_elManual171_AThr7_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_LEDGF_myotube_filtered_unireads_elManual171_AThr5_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_RNAPII_myoblast_filtered_unireads_elEst181_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_RNAPII_myotube_filtered_unireads_elEst171_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_SPT16_myoblast_filtered_unireads_elEst181_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_SPT16_myotube_filtered_unireads_elEst141_AThr3_bin50-RPM_BGSub-scaleReverse-Spikedin.bw",
        "WIGfs_HGDF2_KO_1-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_1-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_1-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_1-D6-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_CageMTresc-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_CageMTresc-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_CageMTresc-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_CageMTresc-D6-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2-D6-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_WTresc-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_WTresc-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_WTresc-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_2_WTresc-D6-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_3-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_3-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_3-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_HGDF2_KO_3-D6-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_WT-D0-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_WT-D0-rep2_unireads_elManual0_bin50.bw",
        "WIGfs_WT-D6-rep1_unireads_elManual0_bin50.bw",
        "WIGfs_WT-D6-rep2_unireads_elManual0_bin50.bw")


expname_vec <- c("H3K27me3_myoblast",
        "H3K27me3_myotube",
        "H3K36me2_myoblast",
        "H3K36me2_myotube",
        "H3K36me3_myoblast",
        "H3K36me3_myotube",
        "HDGF2_myoblast",
        "HDGF2_myotube",
        "LEDGF_myoblast",
        "LEDGF_myotube",
        "RNAPII_myoblast",
        "RNAPII_myotube",
        "SPT16_myoblast",
        "SPT16_myotube",
        "HGDF2_KO_1-D0-rep1",
        "HGDF2_KO_1-D0-rep2",
        "HGDF2_KO_1-D6-rep1",
        "HGDF2_KO_1-D6-rep2",
        "HGDF2_KO_2_CageMTresc-D0-rep1",
        "HGDF2_KO_2_CageMTresc-D0-rep2",
        "HGDF2_KO_2_CageMTresc-D6-rep1",
        "HGDF2_KO_2_CageMTresc-D6-rep2",
        "HGDF2_KO_2-D0-rep1",
        "HGDF2_KO_2-D0-rep2",
        "HGDF2_KO_2-D6-rep1",
        "HGDF2_KO_2-D6-rep2",
        "HGDF2_KO_2_WTresc-D0-rep1",
        "HGDF2_KO_2_WTresc-D0-rep2",
        "HGDF2_KO_2_WTresc-D6-rep1",
        "HGDF2_KO_2_WTresc-D6-rep2",
        "HGDF2_KO_3-D0-rep1",
        "HGDF2_KO_3-D0-rep2",
        "HGDF2_KO_3-D6-rep1",
        "HGDF2_KO_3-D6-rep2",
        "WT-D0-rep1",
        "WT-D0-rep2",
        "WT-D6-rep1",
        "WT-D6-rep2")


folder_on_server <- c("gary/myoblast_myotube");


reference_genome <- c("hg19");

output_file <- "/home/descostes/Documents/analysis/fact_ledgf/myoblast-myotube-phoenix2.txt";

color_RGB <- c("11,158,119",
        "11,158,119",
        "217,95,2",
        "217,95,2",
        "117,112,179",
        "117,112,179",
        "231,41,138",
        "231,41,138",
        "102,166,30",
        "102,166,30",
        "230,171,2",
        "230,171,2",
        "166,118,29",
        "166,118,29",
        rep("102,102,102",24))



##############
# MAIN
##############

if(length(bigwig_file_name_vec) != length(expname_vec))
{
    stop("\n bigwig names and expvec should have the same length");
}


to_write <- vector();

for(i in 1:length(bigwig_file_name_vec)) 
{
    current_bigwig <- bigwig_file_name_vec[i];
    current_expname <- expname_vec[i];
    current_color <- color_RGB[i]
    
    to_write[i] <- paste("track type=bigWig name=\"", current_expname,"\" bigDataUrl=http://www.hpc.med.nyu.edu/~descon01/", folder_on_server, "/", current_bigwig, 
            " visibility=full autoScale=off viewLimits=0:50 alwaysZero=on windowingFunction=mean maxHeightPixels=50:50:50, color=", current_color, " db=", reference_genome, sep="");
}


write(to_write, file=output_file, ncolumns=1);


