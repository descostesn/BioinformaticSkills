
bam_files_vec <- c("/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/dKO_C6_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/dKO_C6_rep2_cardioD5_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/dKO_C6_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/dKO_C6_rep2_N2B27-D8_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/HDGF2KO_C6_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/HDGF2KO_C6_rep2_cardioD5_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/HDGF2KO_C6_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/HDGF2KO_C6_rep2_N2B27-D8_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_F3_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_F3_rep2_cardioD5_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_F3_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_F3_rep2_N2B27-D8_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_G8_rep1_cardioD5_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_G8_rep2_cardioD5_SORTED_PICARD_COOR.bam",
        "/ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/WT_rep2_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_G8_rep1_N2B27-D8_SORTED_PICARD_COOR.bam /ifs/home/descon01/data/data_june2017/bam_files/ledgf_rnaseq/PSIPKO_G8_rep2_N2B27-D8_SORTED_PICARD_COOR.bam");


single_end <- "TRUE"

samples_info_file <- c("/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_dKOC6_cardioD5_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_dKOC6_N2B27_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_HDGF2KOC6_cardioD5_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_HDGF2KOC6_N2B27_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_PSIPKOF3_cardioD5_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_PSIPKOF3_N2B27_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_PSIPKOG8_cardioD5_june2017.csv",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/info_files/june2017/comparisonWT_PSIPKOG8_N2B27_june2017.csv")


name_of_reference <- c("ctrl")

output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_dKOC6_cardioD5_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_dKOC6_N2B27_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_HDGF2KOC6_cardioD5_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_HDGF2KOC6_N2B27_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_PSIPKOF3_cardioD5_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_PSIPKOF3_N2B27_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_PSIPKOG8_cardioD5_june2017/",
        "/ifs/home/descon01/analysis/fact_ledgf/rnaseq/diff_rna_seq/DESeq2/comparisonWT_PSIPKOG8_N2B27_june2017/")


comparison_name <- c("comparisonWT_dKOC6_cardioD5_june2017",
        "comparisonWT_dKOC6_N2B27_june2017",
        "comparisonWT_HDGF2KOC6_cardioD5_june2017",
        "comparisonWT_HDGF2KOC6_N2B27_june2017",
        "comparisonWT_PSIPKOF3_cardioD5_june2017",
        "comparisonWT_PSIPKOF3_N2B27_june2017",
        "comparisonWT_PSIPKOG8_cardioD5_june2017",
        "comparisonWT_PSIPKOG8_N2B27_june2017")


refseq_anno <- "/ifs/home/descon01/cluster/Annotations/mouse/mm10/gtf_files/feb2017/refseq_mm10-newfile-filtered.gtf"


expnames_vec <- c("WT_rep1 WT_rep2 dKO_C6_rep1 dKO_C6_rep2",
        "WT_rep1 WT_rep2 dKO_C6_rep1 dKO_C6_rep2",
        "WT_rep1 WT_rep2 HDGF2KO_C6_rep1 HDGF2KO_C6_rep2",
        "WT_rep1 WT_rep2 HDGF2KO_C6_rep1 HDGF2KO_C6_rep2",
        "WT_rep1 WT_rep2 PSIPKO_F3_rep1 PSIPKO_F3_rep2",
        "WT_rep1 WT_rep2 PSIPKO_F3_rep1 PSIPKO_F3_rep2",
        "WT_rep1 WT_rep2 PSIPKO_G8_rep1 PSIPKO_G8_rep2",
        "WT_rep1 WT_rep2 PSIPKO_G8_rep1 PSIPKO_G8_rep2")

output_format <- "pdf";

sva <- FALSE

to_write <- paste(bam_files_vec, single_end, samples_info_file, name_of_reference, output_folder, comparison_name, refseq_anno, expnames_vec, output_format, sva, sep=";");

write(to_write, file="ledgf_rnaseq_june2017_pdfredo.conf", ncolumns=1);












