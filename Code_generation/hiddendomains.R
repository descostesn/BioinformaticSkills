options(scipen=999)

chromInfo_file <- "/ifs/home/descon01/programs/hiddenDomains/ChromInfo_hg19.txt"
bin_size <- c(50000, 20000, 100000, 10000)

expBam <- c("/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0_HDGF2_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0_PolII_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-SPT16_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5_HDGF2_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5_PolII_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5_SPT16_mm10_ESC.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K27me3_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/H3K27me3_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/HDGF2_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/HDGF2_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/RNAPII_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/RNAPII_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam")


inputbam <- c("/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day0-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/gary_paper_data/GEO_submission/Day5-input_mm10_ESC.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myoblast_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_november2017/bam_files/ledgf_myoblast_myotube/input_myotube_filtered_bwalign_hg19_v3_m1_k1_SORTED_PICARD_COOR.bam")

output_folder_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/ESC/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myoblast/with_hiddenDomains/",
        "/ifs/home/descon01/analysis/fact_ledgf/peak_calling/myotube/with_hiddenDomains/")

expname_vec <- c("Day0_HDGF2_mm10_ESC.bam",
        "Day0_PolII_mm10_ESC.bam",
        "Day0-SPT16_mm10_ESC.bam",
        "Day5_HDGF2_mm10_ESC.bam",
        "Day5_PolII_mm10_ESC.bam",
        "Day5_SPT16_mm10_ESC.bam",
        "H3K27me3_myoblast",
        "H3K27me3_myotube",
        "HDGF2_myoblast",
        "HDGF2_myotube",
        "RNAPII_myoblast",
        "RNAPII_myotube")

to_write <- vector();

for(i in 1:length(bin_size)) 
{
    outputFilePrefix <- paste(output_folder_vec, bin_size[i], "/", expname_vec, sep="");
    
    for(j in 1:length(output_folder_vec)){
        
        outputfoldertocreate <- paste(output_folder_vec[j], bin_size[i], "/", sep="");
        
        if(!file.exists(outputfoldertocreate))
        {
            dir.create(outputfoldertocreate, recursive = TRUE)
        }
    }
    
    to_write <- c(to_write, paste(chromInfo_file, bin_size[i], expBam, inputbam, outputFilePrefix, sep=";"));
}

cat("Length should be ", length(expBam)*length(bin_size))


write(to_write, file="gary_ESCEB_MBMT.conf", ncolumns=1)



