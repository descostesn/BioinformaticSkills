########
# This script makes intersections for an unlimited number of bed files
# Descostes feb2018
#########

library("parallel")
library("GenomicRanges")
options(scipen=999)

################
# PARAMETERS
################

bed_files_path_vec <- c("/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/K27me3_K27M0Hr_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/K27me3_K27WT0Hr_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_12hr-H3K27M-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_12hr-H3WT-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_24hr-H3K27M-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_24hr-H3WT-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_6hr-H3K27M-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_6hr-H3WT-293TREX_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_72H_H33K27M_analysis_100000.bed",
        "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/H3K27me3_72H_H33WT_analysis_100000.bed")

output_folder <- "/ifs/home/descon01/analysis/K27M_project/peak_calling/with_hiddendomains/by_months/feb2018/selection/"

output_filename <- "intersect-H3K27me3.bed"


################



##############
# MAIN
##############

all_extensions <- unique(unlist(lapply(bed_files_path_vec, function(x) strsplit(basename(x), "\\.")[[1]][[2]])))

if(isTRUE(all.equal(length(all_extensions), 1))){
    if(!isTRUE(all.equal(all_extensions,"bed")))
        stop("Only bed format is accepted")
}else stop("All files must be of the same format")


cat("Reading bed files\n");

bed_table_list <- mclapply(bed_files_path_vec, function(x) return(read.table(x, stringsAsFactors=FALSE)))
bed_table <- do.call(rbind,bed_table_list)

cat("performing intersection by chromosome\n")

bed_table_granges <- GRanges(seqnames=bed_table[,1], ranges=IRanges(start=bed_table[,2], end=bed_table[,3]))
rle_coverage_list <- coverage(bed_table_granges)

result_bed_list <- list()
result_bed_list <- mcmapply(function(x, chr_name){
            ## Try this code with: 
            ## gr1 <- GRanges(seqnames=c("chr1", "chr1", "chr1","chr1","chr1"), ranges=IRanges(start=c(1,3,20, 23, 22), end=c(5,6,25, 27, 26)))
            ## x <- coverage(gr1)$chr1
            ## it should remain the intervals 3-5 and 23-25
            value_chr <- as.numeric(runValue(x))
            start_chr <- as.numeric(start(x)[which(value_chr == length(bed_files_path_vec))])
            end_chr <- as.numeric(end(x)[which(value_chr == length(bed_files_path_vec))])
            to_keep <- (end_chr-start_chr) > 0
            
            if(length(which(to_keep)))
                return(data.frame(chr= chr_name, start=start_chr[to_keep], end=end_chr[to_keep], name="Intersect"))
            else return(NA)
        }, rle_coverage_list, unique(bed_table[,1]), SIMPLIFY=FALSE)

result_bed <- na.omit(do.call(rbind, result_bed_list))

cat("Writing intersected file\n")
write.table(result_bed, file=paste(output_folder, output_filename,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


