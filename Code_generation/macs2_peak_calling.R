#BAMFILEVEC;INPUTFILEVEC;EXPNAME;OUTPUTFOLDER;FORMAT;GENOMESIZE;TAGSIZE;PVALUE;ELONGATIONSIZE;ARTEFACTTHRESHOLD


bamfile_vec <- c("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K27me3_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K27me3_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K36me2_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K36me2_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K36me3_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K36me3_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/HDGF2_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/HDGF2_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/LEDGF_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/LEDGF_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/SPT16_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/SPT16_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam");

input_file_vec <- c("/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam");


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
        "SPT16_myoblast",
        "SPT16_myotube");


output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_macs2/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_macs2/");
        

format <- c("BAM");

genome_size <- c(2.7e9);

tagsize_vec <- rep(51,12)

artefact_threshold <- c(5,4,6,5,6,6,5,4,7,5,3,3)

elongation_size <- c(191,191,163,163,115,115,181,181,171,171,181,141)

q_value <- c(0.001,  0.01,  0.02,  0.03,  0.04,  1e-04,  1e-05);


if(length(bamfile_vec) != length(input_file_vec) || length(input_file_vec) != length(expname_vec) || length(expname_vec) != length(tagsize_vec) || length(tagsize_vec) != length(artefact_threshold) || length(artefact_threshold) != length(elongation_size))
{
	stop("\n PROBLEM IN LENGTH OF ELEMENTS")
}

to_write <- vector();
count <- 1;

for(i in 1:length(bamfile_vec))
{
	for(j in 1:length(q_value))
	{
		to_write[count] <- paste(bamfile_vec[i], input_file_vec[i], expname_vec[i], paste(output_folder[i], q_value[j], "/", sep=""), format, genome_size, tagsize_vec[i], q_value[j], elongation_size[i], artefact_threshold[i], sep=";");
		count <- count+1;
	}
}


cat("The number of lines should be: ", length(bamfile_vec)*length(q_value));

write(to_write, file="myoblastmyotube_nov2017_minuspolII.conf", ncolumns=1)






