##############
# This script performs time course DE analysis of RNA-seq data with the bioconductor package maSigPro.
# Descostes May 2017.
##############



library("maSigPro");
library("Rsubread");
library("edgeR");
library("NGSprofiling");
library("Rargs");
library("pheatmap");
library("RColorBrewer");

################
# PARAMETERS
################


paramsDefinition <- list();


# Required arguments
paramsDefinition[["--bamFilesVec"]] <- list(variableName="bam_files_vec", numeric=F, mandatory=T, description="Vector of bam files from which the different clusters will be plotted.");
paramsDefinition[["--refseqAnno"]] <- list(variableName="refseq_anno", numeric=F, mandatory=T, description="GTF files of the ant annotations.");
paramsDefinition[["--expnamesVec"]] <- list(variableName="expnames_vec", numeric=F, mandatory=T, description="Vector of expnames corresponding to the bam files.");
paramsDefinition[["--designFile"]] <- list(variableName="design_file", numeric=F, mandatory=T, description="Path to the design file.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the output folder.");
paramsDefinition[["--minClusterNb"]] <- list(variableName="min_cluster_nb", numeric=T, mandatory=T, description="Numeric giving the minimum number of cluster to start with. Go up to 10.");


#bam_files_vec <- c("/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/W_rep1_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/W_rep2_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/W_rep3_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/W_rep4_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/G_rep2_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/G_rep3_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/G_rep4_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/RG_rep2_fatbodyAligned.sortedByCoord.out.bam",
#		"/ifs/home/descon01/data/data_april2017/bam_files/comzit_ant_fatbody/RG_rep3_fatbodyAligned.sortedByCoord.out.bam");
#
#refseq_anno <- "/ifs/home/descon01/cluster/Annotations/ants/Hsalv3_5/gtf_files/hsal.v3_5.gtf";
#
#expnames_vec <- c("W_rep1",
#		"W_rep2",
#		"W_rep3",
#		"W_rep4",
#		"G_rep2",
#		"G_rep3",
#		"G_rep4",
#		"RG_rep2",
#		"RG_rep3");
#
#design_file <- "/ifs/home/descon01/analysis/comzit_ants/differential_analysis/info_file/fatbody/filtered/timecourse_W_G_RG_1group.txt";
#
#output_folder <- "/ifs/home/descon01/analysis/comzit_ants/differential_analysis/results/maSigPro/fatbody_W_G_RG/";
#
#min_cluster_nb <- 4;

################



################
# FUNCTION
################


find_matrix_per_group <- function(result, edgeR_cpm_matrix){
    
    #Create lists of gene names per group
    groups_list <- split(result$cut, factor(result$cut));
    
    # Retrieve rows in the input matrix per group
    matrix_per_group_list <- lapply(groups_list, function(x, input_matrix){
                
                names_vec <- names(x);
                indexes <- match(names_vec, rownames(input_matrix));
                
                if(length(which(is.na(indexes))) != 0)
                {
                    stop("\n some genes name were not retrieved\n");
                    
                }
                
                if(length(names_vec) != length(indexes))
                {
                    stop("\n You do not have a one to one relation in your structure\n");
                }
                
                return(input_matrix[indexes,]);
                
            }, edgeR_cpm_matrix)
    
    return(matrix_per_group_list);
}


output_results <- function(result, edgeR_cpm_matrix, matrix_output_folder){

    matrix_per_group_list <- find_matrix_per_group(result, edgeR_cpm_matrix);
	
	#Write all matrice
	checkingOutputFolder(matrix_output_folder);
	matrix_nb <- 1;
	write_tables <- lapply(matrix_per_group_list, function(x, out){
				write.table(x, file=paste0(out, "cluster_", matrix_nb, ".txt"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
				matrix_nb <<- matrix_nb + 1;
				return(NULL);
			}, matrix_output_folder);
}


boxplot_per_group <- function(result, edgeR_cpm_matrix, output_folder){
    
    matrix_per_group_list <- find_matrix_per_group(result, edgeR_cpm_matrix);
    output_folder_boxplots <- paste0(output_folder, "boxplots/");
    checkingOutputFolder(output_folder_boxplots);
    
    matrix_nb <- 1;
    
    lapply(matrix_per_group_list, function(mat,out){
                
                if(!is.null(ncol(mat))) # otherwise matrix only contains one line
                {
                    colors_tab <- brewer.pal(ncol(mat),"Set3");
                    
                    png(filename=paste(out, "cluster", matrix_nb, "-notches.png",sep=""), width = 600, height = 600, bg = "transparent")
                    par(mfrow=c(1,2))
                    boxplot(mat, ylab= "TPM values", col= colors_tab, outline=FALSE, names = colnames(mat), las=2, cex.axis=0.7, notch=TRUE);
                    boxplot(mat, ylab= "TPM values", col= colors_tab, outline=TRUE, names = colnames(mat), las=2, cex.axis=0.7, notch=TRUE);
                    dev.off();
                    
                    png(filename=paste(out, "cluster", matrix_nb, ".png",sep=""), width = 600, height = 600, bg = "transparent")
                    par(mfrow=c(1,2))
                    boxplot(mat, ylab= "TPM values", col= colors_tab, outline=FALSE, names = colnames(mat), las=2, cex.axis=0.7, notch=FALSE);
                    boxplot(mat, ylab= "TPM values", col= colors_tab, outline=TRUE, names = colnames(mat), las=2, cex.axis=0.7, notch=FALSE);
                    dev.off();
                    
                }
                
                matrix_nb <<- matrix_nb + 1;
                
            }, output_folder_boxplots)
    
}

################





##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


#control
		
if(length(bam_files_vec) != length(expnames_vec))
{
	stop("One experiment name should be given per bam file\n");
}

checkingOutputFolder(output_folder);
log_report <- vector();

# Building the count table from experiments

cat("Computing the input counts matrix\n");
result_counts <- featureCounts(bam_files_vec, annot.ext= refseq_anno, isGTFAnnotationFile=TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", useMetaFeatures=TRUE, allowMultiOverlap=FALSE, isPairedEnd=FALSE);	
counts_matrix <- result_counts$counts;
colnames(counts_matrix) <- expnames_vec;

# Removing the lines having too many zero (this can be a problem for the fit)

cat("Keeping only genes having more than 2 reads on average\n");
mean_nb_reads <- as.numeric(unlist(apply(counts_matrix, MARGIN=1,mean)));
to_keep <- which(mean_nb_reads > 2);
percent_kept <- (length(to_keep)*100)/nrow(counts_matrix);
cat("\t Keeping ", length(to_keep), "/", nrow(counts_matrix), "(", percent_kept, "%)");

counts_matrix <- counts_matrix[to_keep,];



# Normalization by TPM

cat("TPM normalization\n");
edgeR_cpm_matrix <- cpm(counts_matrix);


# Reading the design file

edesign_df <- read.table(design_file, header=TRUE, sep=";", row.names = 1);


# Creation of the regression matrix and making a regression fit for each gene.

cat("Making the regression matrix\n");

design <- make.design.matrix(edesign_df, degree = length(unique(edesign_df$Time))-1);

cat("Making a regression fit for each gene\n");

fit <- p.vector(edgeR_cpm_matrix, design, Q = 0.05, MT.adjust = "BH", counts=TRUE, family= negative.binomial(10));

# Recording statistics on the fit
log_report <- c(log_report, "Statistics on the fit:");
log_report <- c(log_report, paste0("\t Number of input genes: ", fit$G));
log_report <- c(log_report, paste0("\t Number of significantly varying genes: ", fit$i));
log_report <- c(log_report, paste0("\t FDR: ", fit$FDR));


cat("Finding significant differences and retrieving the genes\n");

tstep <- T.fit(fit, step.method = "backward", alfa = 0.05, family = negative.binomial(10));


cat("Plotting different summaries of the analysis\n");

sigs <- get.siggenes(tstep, rsq = 0.7, vars = "groups");

png(filename=paste(output_folder, "venndiagram.png",sep=""), width = 600, height = 600, bg = "transparent");
suma2Venn(sigs$summary);
dev.off();


# Generating clusters		
		
sigs <- get.siggenes(tstep, rsq = 0.7, vars = "all");

write.table(cbind(sigs$sig.genes$sig.profiles, sigs$sig.genes$sig.pvalues), file=paste0(output_folder, "siggenes_table.txt"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);

#test <- see.genes(sigs$sig.genes, show.fit=T, dis=design$dis, cluster.method="hclust", cluster.data=1, k=3)

count_computation <- 1;
cluster_data_names <- c("sig.profiles","coefficients","t.score");
total_nb_operation_hclust <- 3*length(min_cluster_nb:10)*3*3;
total_nb_operation_kmean <- 3*length(min_cluster_nb:10);
total_nb_operation <- total_nb_operation_hclust + total_nb_operation_kmean;

for(data_type in 1:2)
{
	for(nb_cluster in min_cluster_nb:10)
	{
		for(cluster_method in c("hclust", "kmeans"))
		{
				for(summary_mode in c("median")) #c("representative", "median"), remove representative because of a bug. Posted on bioconductor but no reply so far
				{
					if(cluster_method == "hclust")
					{
						for(agglo_method in c("single", "complete", "average"))
						{
							for(distance_metrics in c("cor", "euclidean", "manhattan"))
							{
								cat("Computing ", count_computation, "/", total_nb_operation,"\n");
								
																
								pdf(file=paste0(output_folder, cluster_data_names[data_type], "_", 
												nb_cluster, "-clusters_",
												cluster_method,
												summary_mode, "_",
												agglo_method, "_",
												distance_metrics, ".pdf"), width=10, height=10, onefile=TRUE)
								result <- see.genes(data=sigs$sig.genes, 
										edesign = design$edesign, 
										time.col = 1, 
										repl.col = 2, 
										group.cols = c(3:ncol(design$edesign)), 
										names.groups = colnames(design$edesign)[3:ncol(design$edesign)],
										cluster.data = data_type, 
										groups.vector = design$groups.vector, 
										k = nb_cluster, 
										cluster.method = cluster_method, 
										distance = distance_metrics,
										agglo.method = agglo_method,
										show.fit = FALSE,
										dis = design, 
										step.method = "backward", 
										min.obs = 3,
										alfa = 0.05,  
										show.lines = TRUE,
										iter.max = 500,
										summary.mode = summary_mode, 
										newX11 = FALSE);
								dev.off();
								count_computation <<- count_computation + 1;

								matrix_output_folder <- paste0(output_folder, cluster_data_names[data_type], "_", nb_cluster, "-clusters_", cluster_method, summary_mode, "_", agglo_method, "_", distance_metrics, "/");
								output_results(result, edgeR_cpm_matrix, matrix_output_folder);
        boxplot_per_group(result, edgeR_cpm_matrix, matrix_output_folder);
        
								#Generating a heatmap of the matrix
								
							   pheatmap(mat = sigs$sig.genes[[cluster_data_names[data_type]]],
									   scale= "row",
									   cluster_cols = FALSE,
									   clustering_distance_rows = if(distance_metrics == "cor") "correlation" else distance_metrics,
									   clustering_method = agglo_method,
									   cutree_rows = nb_cluster, 
									   annotation_names_row = FALSE, 
									   annotation_names_col = TRUE,
									   drop_levels = FALSE, 
									   show_rownames = FALSE,  
									   filename = paste0(output_folder, cluster_data_names[data_type], "_", nb_cluster, "-clusters_", cluster_method, summary_mode, "_", agglo_method, "_", distance_metrics, "heatmap.pdf"), 
									   silent = TRUE);
							}
						}
						
					}else{
						cat("Computing ", count_computation, "/", total_nb_operation,"\n");
												
						pdf(file=paste0(output_folder, cluster_data_names[data_type], "_", 
										nb_cluster, "-clusters_",
										cluster_method,
										summary_mode, ".pdf"), width=10, height=10, onefile=TRUE)
	
						result <- see.genes(data=sigs$sig.genes, 
								edesign = design$edesign, 
								time.col = 1, 
								repl.col = 2, 
								group.cols = c(3:ncol(design$edesign)), 
								names.groups = colnames(design$edesign)[3:ncol(design$edesign)],
								cluster.data = data_type, 
								groups.vector = design$groups.vector, 
								k = nb_cluster, 
								cluster.method = cluster_method,
								show.fit = FALSE,
								dis = design, 
								step.method = "backward", 
								min.obs = 3, 
								alfa = 0.05,  
								show.lines = TRUE,
								iter.max = 500,
								summary.mode = summary_mode, 
								newX11 = FALSE);
						
						dev.off();
						
						count_computation <<- count_computation + 1;
						
						matrix_output_folder <- paste0(output_folder, cluster_data_names[data_type], "_", nb_cluster, "-clusters_", cluster_method, summary_mode, "/");
						output_results(result, edgeR_cpm_matrix, matrix_output_folder);
					}
					
				}	
			}
		} 
	}









