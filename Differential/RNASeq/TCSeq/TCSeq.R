#########
# This script aims at analyzing time-course rna-seq data with TCSeq.
# Descostes ap 2018
#########

library(TCseq)

#########################
# PARAMETERS
#########################


#bam_files_vec <- c("ovary_Worker_rep1Aligned.sortedByCoord.out.bam",
#        "ovary_Worker_rep2Aligned.sortedByCoord.out.bam",
#        "ovary_Worker_rep3Aligned.sortedByCoord.out.bam",
#        "ovary_Worker_rep4Aligned.sortedByCoord.out.bam",
#        "ovary_Gamergate_rep1Aligned.sortedByCoord.out.bam",
#        "ovary_Gamergate_rep2Aligned.sortedByCoord.out.bam",
#        "ovary_Gamergate_rep3Aligned.sortedByCoord.out.bam",
#        "ovary_Gamergate_rep4Aligned.sortedByCoord.out.bam",
#        "ovary_Reverted_gamergate_rep1Aligned.sortedByCoord.out.bam",
#        "ovary_Reverted_gamergate_rep2Aligned.sortedByCoord.out.bam",
#        "ovary_Reverted_gamergate_rep3Aligned.sortedByCoord.out.bam",
#        "ovary_Reverted_gamergate_rep4Aligned.sortedByCoord.out.bam")


#bam_files_vec <- c("W_rep1_fatbodyAligned.sortedByCoord.out.bam",
#"W_rep2_fatbodyAligned.sortedByCoord.out.bam",
#"W_rep3_fatbodyAligned.sortedByCoord.out.bam",
#"W_rep4_fatbodyAligned.sortedByCoord.out.bam",
#"G_rep1_fatbodyAligned.sortedByCoord.out.bam",
#"G_rep2_fatbodyAligned.sortedByCoord.out.bam",
#"G_rep3_fatbodyAligned.sortedByCoord.out.bam",
#"G_rep4_fatbodyAligned.sortedByCoord.out.bam",
#"RG_rep1_fatbodyAligned.sortedByCoord.out.bam",
#"RG_rep2_fatbodyAligned.sortedByCoord.out.bam",
#"RG_rep3_fatbodyAligned.sortedByCoord.out.bam",
#"RG_rep4_fatbodyAligned.sortedByCoord.out.bam")

bam_files_vec <- c("CentralBrain_Worker_rep1Aligned.sortedByCoord.out.bam",
        "Centrarain_Worker_rep2Aligned.sortedByCoord.out.bam",
        "Centrarain_Worker_rep3Aligned.sortedByCoord.out.bam",
        "Centrarain_Worker_rep4Aligned.sortedByCoord.out.bam",
        "Centrarain_Gamergate_rep1Aligned.sortedByCoord.out.bam",
        "Centrarain_Gamergate_rep2Aligned.sortedByCoord.out.bam",
        "Centrarain_Gamergate_rep3Aligned.sortedByCoord.out.bam",
        "Centrarain_Gamergate_rep4Aligned.sortedByCoord.out.bam",
        "Centrarain_Reverted_gamergate_rep1Aligned.sortedByCoord.out.bam",
        "Centrarain_Reverted_gamergate_rep2Aligned.sortedByCoord.out.bam",
        "Centrarain_Reverted_gamergate_rep3Aligned.sortedByCoord.out.bam",
        "Centrarain_Reverted_gamergate_rep4Aligned.sortedByCoord.out.bam")

bam_dir <- "/home/descostes/Documents/analysis/comzit_ants/bamfiles/"

refseq_anno_path <- "/home/descostes/Documents/analysis/comzit_ants/new_version/hsal_v8.5.gff"

sample_id <- c("W_rep1",
		"W_rep2",
		"W_rep3",
		"W_rep4",
        "G_rep1",
		"G_rep2",
		"G_rep3",
		"G_rep4",
        "RG_rep1",
		"RG_rep2",
		"RG_rep3",
        "RG_rep4")

group <- timepoint <- c("W",
        "W",
        "W",
        "W",
        "G",
        "G",
        "G",
        "G",
        "RG",
        "RG",
        "RG",
        "RG")

minimum_number_of_tags <- 5
minimum_number_exp_with_minNumTags <- 4

# a character vector, each charcter string in the vector gives a contrast of 
# two groups with the format group2vsgroup1', group1 is the denominator level 
# in the fold changes and group2 is the numerator level in the fold changes.
contrasts_vec <- c("GvsW", "RGvsW")

#output_folder <- "/home/descostes/Documents/analysis/comzit_ants/differential_analysis/results/time_course_experiment/second_analysis_ap2018/TCseq/ovary/ztransform/"
#output_folder <- "/home/descostes/Documents/analysis/comzit_ants/differential_analysis/results/time_course_experiment/second_analysis_ap2018/TCseq/fatbody/ztransform/"
output_folder <- "/home/descostes/Documents/analysis/comzit_ants/differential_analysis/results/time_course_experiment/second_analysis_ap2018/TCseq/brain/ztransform/"

z_transform <- TRUE

#########################
# MAIN
#########################


if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)

group_nb_vec <- c(W = 1, G = 2, RG = 3)
group_nb_vec <- group_nb_vec[timepoint]

cat("Converting GFF to interval format\n");

refseq_anno <- read.table(refseq_anno_path, stringsAsFactors=FALSE)
genomic_intervals <- data.frame(id= refseq_anno[,2], chr=refseq_anno[,1], 
        start=refseq_anno[,4], end = refseq_anno[,5])

cat("Building bam data.frame and object\n")

experiment_bamfile <- data.frame(sampleid = sample_id, timepoint = timepoint, 
        group = group, BAMfile = bam_files_vec)
tca <- TCA(design = experiment_bamfile, genomicFeature = genomic_intervals)
tca <- countReads(tca, dir = bam_dir)


cat("Performing the differential analysis\n")

tca <- DBanalysis(tca, filter.type = "raw", 
        filter.value = minimum_number_of_tags, 
        samplePassfilter = minimum_number_exp_with_minNumTags)


cat("\t Two by two\n")

DBres <- DBresult(tca, 
        contrasts = contrasts_vec,
        p.adjust = "BH", 
        top.sig = TRUE, 
        pvalue = "paj",
        pvalue.threshold = 0.05, 
        abs.fold = 2, 
        direction = "both",
        result.type = "list")

mock <- mapply(function(resTable, nameTable, output_fold){
            write.table(resTable, file = paste0(output_fold, nameTable, 
                            "DE_table.txt"), sep="\t", quote=FALSE, 
                    row.names = TRUE, col.names = TRUE)
        }, DBres, contrasts_vec, MoreArgs = list(output_folder))


cat("\t Time course\n")

tca <- timecourseTable(tca, 
        value = "expression", 
        lib.norm = TRUE,
        norm.method = "cpm", filter = TRUE, pvalue = "BH",
        pvalue.threshold = 0.05, abs.fold = 2, direction = "both")

cpm_matrix_counts <- counts(tca, normalization="cpm")
colnames(cpm_matrix_counts) <- sample_id
write.table(cpm_matrix_counts, file = paste0(output_folder, 
                "total_cpm_counts.txt"), sep="\t", quote=FALSE, 
        row.names = TRUE, col.names = TRUE)

cat("\t\t Performing clustering\n")

for(algo_type in c("km", "hc", "cm")){
    
    for(distance_method in c("euclidean", "manhattan")){
        
        for(nb_groups in seq_len(10))
        {
            if(!isTRUE(all.equal(nb_groups,1)) && 
                    !(isTRUE(all.equal(algo_type, "km")) && isTRUE(all.equal(distance_method, "manhattan")))){
                
                cat("\t\t\t Processing: ", algo_type, " - ", distance_method, " - ", nb_groups,"\n")
                
                current_clust <- timeclust(tca, 
                        algo = algo_type, 
                        k = nb_groups, 
                        dist = distance_method,
                        standardize = z_transform)
                
                cluster_output_fold <- paste0(output_folder, 
                        paste(algo_type, distance_method, nb_groups, sep="_"), "/")
                if(!file.exists(cluster_output_fold))
                    dir.create(cluster_output_fold, recursive=TRUE)
                
                png(filename = paste0(cluster_output_fold, "profiles.png"), width = 1000, height = 1000)
                p <- timeclustplot(object = current_clust, 
                        categories = "timepoint",
                        value = paste0("zscore(", algo_type, "_", distance_method),
                        cols = 3)
                dev.off()
                
                for (i in seq_len(length(p))) {
                    
                    result <- p[[i]]
                    output_table <- do.call(rbind,lapply(split(result[[1]],
                                            factor(result[[1]][,1])), 
                                    function(x){return(
                                                paste(x$value, 
                                                        collapse=" "))}))
                    colnames(output_table) <- paste(unique(group), 
                            collapse=" ")
                    write.table(output_table, 
                            file = paste0(cluster_output_fold, "cluster-", i, 
                                    ".txt"), sep="\t", quote=FALSE, 
                            row.names = TRUE, col.names = TRUE)
                    
                    # Retrieving the cpm counts to plot individual genes in cluster folder
                    names_vec <- rownames(output_table)
                    cpm_counts_tmp_mat <- cpm_matrix_counts[
                            sapply(names_vec, function(x){
                                        which(rownames(cpm_matrix_counts) 
                                        == x)}),]
    
                    cluster_folder <- paste0(cluster_output_fold, "cluster-", i, "/")
                    if(!file.exists(cluster_folder))
                        dir.create(cluster_folder, recursive = TRUE)
                    
                    for (j in seq_len(nrow(cpm_counts_tmp_mat))) {	
                        mean_vec <- unlist(lapply(split(cpm_counts_tmp_mat[j,], factor(timepoint, levels = unique(timepoint))),mean))
                        
                        png(filename = paste0(cluster_folder, names_vec[j], ".png"))
                        plot(group_nb_vec, cpm_counts_tmp_mat[j,])
                        lines(mean_vec, col="red")
                        dev.off()
                    }
                }
            }
        }
    }
}



