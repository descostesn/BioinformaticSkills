###############
# This script performs differential analysis using edgeR package. 
# Descostes october 2016
###############


library("edgeR");
library("GenomicAlignments");
library("GenomicFeatures");
library("statmod");
library("NGSprofiling");
library("Rsubread");
library("Rargs");
library("ggplot2");
library("pheatmap");





################
# PARAMETERS
################



# Define the required options
paramsDefinition=list();

paramsDefinition[["--bamFilesVec"]]=list(variableName="bam_files_vec", numeric=F, mandatory=T, description="A space separated list of bam files to compare.");
paramsDefinition[["--refseqAnno"]]=list(variableName="refseq_anno", numeric=F, mandatory=T, description="File path of the file containing annotations. This should be in gtf format");
paramsDefinition[["--singleEnd"]]=list(variableName="single_end", numeric=F, mandatory=T, description="Boolean indicating if the data are single or paired ended.", postConversion=as.logical);
paramsDefinition[["--expnamesVec"]]=list(variableName="expname_vec", numeric=F, mandatory=T, description="A space separated vector of exp names corresponding to the bam files.");
paramsDefinition[["--conditionsVec"]]=list(variableName="conditions_vec", numeric=T, mandatory=T, description="Vector of numbers indicating the group to which each bam files belong.");
paramsDefinition[["--outputFolder"]]=list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the output folder.");
paramsDefinition[["--species"]]=list(variableName="species", numeric=F, mandatory=T, description="Single string indicating species. should be human, mouse, drosophila or ant");
paramsDefinition[["--gffFile"]]=list(variableName="gff_file", numeric=F, mandatory=T, description="Single string indicating gff file corresponding to the gtf file.");
paramsDefinition[["--groupnames"]]=list(variableName="groupnames_vec", numeric=F, mandatory=T, description="String vector indicating the category names of the bam files.");
paramsDefinition[["--cell"]]=list(variableName="cell_vec", numeric=F, mandatory=T, description="String vector indicating the cell for each bam file.");



#bam_files_vec <- c("/home/descostes/Documents/analysis/comzit_ants/bamfiles/CentralBrain_Worker_rep1Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Worker_rep2Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Worker_rep3Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Worker_rep4Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Gamergate_rep1Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Gamergate_rep2Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Gamergate_rep3Aligned.sortedByCoord.out.bam",
#        "/home/descostes/Documents/analysis/comzit_ants/bamfiles/Centrarain_Gamergate_rep4Aligned.sortedByCoord.out.bam")
#refseq_anno <- "/home/descostes/Documents/analysis/comzit_ants/new_version/hsal_v8.5.gtf"
#single_end <- T
#expname_vec <- c("Wrep1","Wrep2","Wrep3","Wrep4","Grep1","Grep2","Grep3","Grep4")
#conditions_vec <- c(1, 1, 1, 1, 2, 2, 2, 2)
#output_folder <- "/home/descostes/Documents/test/"
#species <- "ant"
#gff_file <- "/home/descostes/Documents/analysis/comzit_ants/new_version/hsal_v8.5.gff"




################



################
# FUNCTION
################


perform_pca <- function(matrix_count, matrix_name, output_folder, groupnames)
{
    pca_result <- prcomp(t(matrix_count));
    
    png(filename=paste(output_folder, "variance_vs_component", matrix_name, ".png",sep=""), width = 600, height = 600, bg = "transparent")
    plot(pca_result, type = "l")
    dev.off();
    
    df_for_plot <- data.frame(pca_result$x, group= groupnames, name=rownames(t(matrix_count)));
    percentVar <- pca_result$sdev^2 / sum( pca_result$sdev^2 );
    
    ggplot(data=df_for_plot, aes_string(x= "PC1", y= "PC2", color="group", label= "name")) + geom_point(size=3) +
            geom_text(aes(y = PC2 + 0.5), position = position_dodge(0.9), vjust = 0.3, hjust=0.7, angle=45) +
            xlab(paste0("PC1", ": ",round(percentVar[1] * 100),"% variance")) +
            ylab(paste0("PC2", ": ",round(percentVar[2] * 100),"% variance")) +
            coord_fixed();
    
    ggsave(paste(matrix_name, "_pca.png",sep=""), plot = last_plot(), device = "png", path = output_folder);
}



plot_fold_change <- function(input_table, output_folder, comparison_title)
{
    fold_change_down <- 2^(-1*input_table$table$logFC[which(input_table$table$logFC<0)]);
    fold_change_up <- 2^(input_table$table$logFC[which(input_table$table$logFC>=0)]);
    
    png(filename=paste(output_folder, comparison_title,".png",sep=""), width = 600, height = 600, bg = "transparent")
    hist(c(fold_change_down, fold_change_up), breaks=100, main="Fold change", xlab="fold change");
    dev.off();
    
}



perform_hiearchical <- function(input_matrix, isDE_bool, cell_vec, groupnames_vec, outputfolder, title_test, title_norm)
{
    de_matrix_counts <- input_matrix[isDE_bool,];
    
    df <- data.frame(cell_vec, groupnames_vec);
    rownames(df) <- colnames(de_matrix_counts);
    output_hierarchical <- paste(outputfolder, "hierarchical/", title_test, "/", title_norm, "/", sep="");
    checkingOutputFolder(output_hierarchical);
    
    nb_row_noVariablity <- which(as.numeric(apply(de_matrix_counts, 1, var))==0);
    
    #if the number of rows without variability is not null, these rows are removed as they are not informative for the heatmap
    
    if(length(nb_row_noVariablity) != 0)
    {
        de_matrix_counts <- de_matrix_counts[-nb_row_noVariablity,];
    }
    
    if(nrow(de_matrix_counts) > 2)
    {
        for(method_clustering in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
        {
            comparison_title <- paste(title_test, title_norm, method_clustering, sep="_");
            png(filename=paste(output_hierarchical, comparison_title,".png",sep=""), width = 600, height = 600, bg = "transparent")
            result_clustering <- pheatmap(de_matrix_counts, annotation=df, clustering_method=method_clustering, cluster_cols = FALSE, treeheight_row=100, show_rownames= FALSE, main= comparison_title, scale= "row");
            dev.off();
            
            write(result_clustering$tree_row$order, file= paste(output_hierarchical, comparison_title, "-order.txt", sep=""), ncolumns=1);
        }
    }
}


return_table_up_or_down <- function(boolean_up_or_down, number_up_or_down, result_test)
{
    if(length(which(unique(boolean_up_or_down) == number_up_or_down)) != 0)
    {
        result_table <- result_test$table[which(boolean_up_or_down == number_up_or_down),];
        result_table <- result_table[order(result_table$logFC, decreasing = if(number_up_or_down > 0) TRUE else FALSE),];
        
    }else{
        result_table <- NA;
    }
    return(result_table);
}


################



##############
# MAIN
##############

# Retreives the parameters
getParams(paramsDefinition);

checkingOutputFolder(output_folder);

if(species != "human" && species != "mouse" && species != "drosophila" && species != "ant")
{
    stop("Species should be human, mouse, drosophila or ant\n\n");
}

###########
# Part1: Processing data and performing differential expression
###########



if(length(unique(conditions_vec)) != 2)
{
    stop("\n This script handles two conditions only\n");
}


if(length(conditions_vec) != length(bam_files_vec))
{
    stop("\n A group nb should be indicated for each bam file\n");
}

if(species != "ant")
{
    cat("Creating a BamFileList\n");
    bamfiles <- BamFileList(bam_files_vec);
    
    cat("Creating the genes database from biomarRt using ensembl\n");
    txdb <- makeTxDbFromGFF(refseq_anno, format="gtf");
    
    cat("Making a Granges list of all exons by genes\n");
    ebg <- exonsBy(txdb, by="gene");
    
    cat("Computing the read counts\n\n");
    se <- summarizeOverlaps(features=ebg, reads=bamfiles,
            mode="Union",
            singleEnd=single_end,
            ignore.strand= TRUE,
            fragments= if(single_end) FALSE else TRUE);
    colnames(se) <- expname_vec;
    counts_matrix <- assay(se);
    
    cat("Preparing the data (conversion, filtering and normalization\n");
    y <- DGEList(counts= counts_matrix, group=conditions_vec);
    
}else{
    
    cat("Computing the input counts matrix\n");
    result_counts <- featureCounts(bam_files_vec, annot.ext= refseq_anno, isGTFAnnotationFile=TRUE, GTF.featureType="CDS", GTF.attrType="gene_id", useMetaFeatures=TRUE, allowMultiOverlap=FALSE, isPairedEnd=FALSE)	
    counts_matrix <- result_counts$counts;
    colnames(counts_matrix) <- expname_vec;
    
    cat("Preparing the data (conversion, filtering and normalization\n");
    y <- DGEList(counts= counts_matrix, group=conditions_vec);
}

#distances on the plot approximate the typical log2 fold changes between the samples
colors_vec <- factor(groupnames_vec);
levels(colors_vec) <- rainbow(length(unique(groupnames_vec)));
png(filename=paste(output_folder, "log2FCplot.png",sep=""), width = 600, height = 600, bg = "transparent");
plotMDS(y, col=as.vector(colors_vec), labels=expname_vec);
dev.off();



keep <- rowSums(cpm(y)>1) >= 2; #keeping genes having at leat 6-7 reads. cpm= count per million
y <- y[keep, , keep.lib.sizes=FALSE];
y <- calcNormFactors(y);
y_classic <- estimateDisp(y);

#Plot the biological coefficient of variation for yclassic
png(filename=paste(output_folder, "biologicalCoefficientOfVariation_exactTest.png",sep=""), width = 600, height = 600, bg = "transparent")
plotBCV(y_classic);
dev.off();


design <- model.matrix(~conditions_vec);
y_glm <- estimateDisp(y, design);

#Plot the biological coefficient of variation for y_glm
png(filename=paste(output_folder, "biologicalCoefficientOfVariation_lrt_qlf_Test.png",sep=""), width = 600, height = 600, bg = "transparent")
plotBCV(y_glm);
dev.off();

cat("\n Performing DE analysis with exact test\n");
exactest_results <- exactTest(y_classic);

cat("\n Performing DE analysis with glm likelihood ratio test\n");
fit_glm <- glmFit(y_glm, design);
results_lrt <- glmLRT(fit_glm, coef=2);

cat("\n Performing DE analysis with glm  QL F-test\n");
fit_qlf <- glmQLFit(y_glm,design)
results_qlf <- glmQLFTest(fit_qlf,coef=2)

#Plot the biological coefficient of variation in the quasi-likelihood case
png(filename=paste(output_folder, "biologicalCoefficientOfVariation_qlf_Test.png",sep=""), width = 600, height = 600, bg = "transparent")
plotQLDisp(fit_qlf);
dev.off();

#Performing PCA on RPKM matrix and not-transformed matrix

cat("Performing PCA\n");

if(species == "ant")
{
    #Retrieving length
    gff_table <- read.table(gff_file, stringsAsFactors=F);
    
    if(nrow(gff_table) != nrow(counts_matrix))
    {
        stop("Problem retrieving gene length\n");
    }
    
    id_to_check <- gff_table$V2;
    index_ordering_gff <- match(rownames(counts_matrix), id_to_check);
    
    if(is.unsorted(index_ordering_gff))
    {
        if(length(which(is.na(match(id_to_check, rownames(counts_matrix))))) != 0)
            stop("\n Problem retrieving gene lengths 2\n");
        
        gff_table <- gff_table[index_ordering_gff,];
    }
    
    length_vec <- gff_table$V5 - gff_table$V4;
    
    
}else{
    
    #Retrieving length
    gff_table <- read.table(gff_file, stringsAsFactors=F);
    
    if(length(which(is.na(match(rownames(counts_matrix), gff_table$V3)))))
    {
        stop("Names of count matrix not found in gff file, did you use concording gtf and gff?\n\n");
    }
    
    index_sorting_gff <- match(rownames(counts_matrix), gff_table$V3);
    gff_table <- gff_table[index_sorting_gff,];
    
    length_vec <- gff_table$V5 - gff_table$V4;
}

counts_matrix_logCPM <- cpm(y, prior.count=2, log=TRUE);
counts_matrix_rpkm <- rpkm(counts_matrix, gene.length= length_vec, log=TRUE);

write.table(counts_matrix, file=paste(output_folder, "count.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F);
write.table(counts_matrix_rpkm, file=paste(output_folder, "count_logrpkm.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F);
write.table(counts_matrix_logCPM, file=paste(output_folder, "count_logCPM.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F);

perform_pca(counts_matrix, "count_matrix", output_folder, groupnames_vec);
perform_pca(counts_matrix_rpkm, "count_matrix_logrpkm", output_folder, groupnames_vec);
perform_pca(counts_matrix_logCPM, "count_matrix_logCpm", output_folder, groupnames_vec);


cat("plotting read counts distribution\n");

png(filename=paste(output_folder, "counts_distribution.png",sep=""), width = 600, height = 600, bg = "transparent")
hist(counts_matrix, breaks=10000, xlim=c(0,1000), xlab="nb_reads (lim=1000)");
dev.off();

cat("Plotting fold change\n");

plot_fold_change(exactest_results, output_folder, "exacTestFoldChange");
plot_fold_change(results_lrt, output_folder, "lrtFoldChange");
plot_fold_change(results_qlf, output_folder, "qlfFoldChange");

cat("Calculating FDR\n");
FDR_exact <- p.adjust(exactest_results$table$PValue, method="BH");
FDR_lrt <- p.adjust(results_lrt$table$PValue, method="BH");
FDR_qlf <- p.adjust(results_qlf$table$PValue, method="BH");

#Including FDR in tables
exactest_results$table <- cbind(exactest_results$table, FDR=FDR_exact);
results_lrt$table <- cbind(results_lrt$table, FDR=FDR_lrt);
results_qlf$table <- cbind(results_qlf$table, FDR=FDR_qlf);

cat("Plotting the libraries size\n");
png(filename=paste(output_folder, "librariessize.png",sep=""), width = 600, height = 600, bg = "transparent");
barplot(y$samples$lib.size*1e-6, names=rownames(y$samples), ylab="Library size (millions)", las=2);
dev.off();


###########
# Part2: Exploring results
###########


cat("\n\n Computing barplots of nb of DE genes for exact test\n");


exact_nb_signif_twofolds <- decideTestsDGE(exactest_results, p.value=0.05, lfc=1)
results_lrt_nb_signif_twofolds <- decideTestsDGE(results_lrt, p.value=0.05, lfc=1)
results_qlf_nb_signif_twofolds <- decideTestsDGE(results_qlf, p.value=0.05, lfc=1)

png(filename=paste(output_folder, "nb_genes_twofolds_pval005_exact.png",sep=""), width = 600, height = 600, bg = "transparent");
barplot(c(summary(exact_nb_signif_twofolds)[3,1], summary(exact_nb_signif_twofolds)[1,1]), names.arg=c(paste("upregulated: ", summary(exact_nb_signif_twofolds)[3,1], sep=""), paste("downregulated: ", summary(exact_nb_signif_twofolds)[1,1],sep="")), main="exact test");
dev.off();

png(filename=paste(output_folder, "nb_genes_twofolds_pval005_lrt.png",sep=""), width = 600, height = 600, bg = "transparent");
barplot(c(summary(results_lrt_nb_signif_twofolds)[3,1],summary(results_lrt_nb_signif_twofolds)[1,1]), names.arg=c(paste("upregulated: ", summary(results_lrt_nb_signif_twofolds)[3,1], sep=""), paste("downregulated: ", summary(results_lrt_nb_signif_twofolds)[1,1],sep="")), main="glm lrt");
dev.off();

png(filename=paste(output_folder, "nb_genes_twofolds_pval005_qlf.png",sep=""), width = 600, height = 600, bg = "transparent");
barplot(c(summary(results_qlf_nb_signif_twofolds)[3,1],summary(results_qlf_nb_signif_twofolds)[1,1]), names.arg=c(paste("upregulated: ", summary(results_qlf_nb_signif_twofolds)[3,1], sep=""), paste("downregulated: ", summary(results_qlf_nb_signif_twofolds)[1,1],sep="")), main="glm qlf");
dev.off();

cat("Generating MA plots\n");

isDE_exact <- as.logical(exact_nb_signif_twofolds);
isDE_lrt <- as.logical(results_lrt_nb_signif_twofolds);
isDE_qlf <- as.logical(results_qlf_nb_signif_twofolds);

DEnames_exact <- rownames(y_classic)[isDE_exact];
DEnames_lrt <- rownames(y_glm)[isDE_lrt];
DEnames_qlf <- rownames(y_glm)[isDE_qlf];

png(filename=paste(output_folder, "exactTest_twofolds_pval005.png",sep=""), width = 600, height = 600, bg = "transparent");
plotSmear(exactest_results, de.tags=DEnames_exact);
abline(h=c(-1,1), col="blue")
dev.off();

png(filename=paste(output_folder, "lrt_twofolds_pval005.png",sep=""), width = 600, height = 600, bg = "transparent");
plotSmear(results_lrt, de.tags=DEnames_lrt);
abline(h=c(-1,1), col="blue")
dev.off();


png(filename=paste(output_folder, "qlf_twofolds_pval005.png",sep=""), width = 600, height = 600, bg = "transparent");
plotSmear(results_qlf, de.tags=DEnames_qlf);
abline(h=c(-1,1), col="blue")
dev.off();

cat("Performing hiearchical clustering\n");

if(length(which(isDE_exact)) > 2)
{
    perform_hiearchical(counts_matrix_logCPM, isDE_exact, cell_vec, groupnames_vec, output_folder, "exact", "CPM_log");
    perform_hiearchical(counts_matrix_rpkm, isDE_exact, cell_vec, groupnames_vec, output_folder, "exact", "RPKM_log");
}

if(length(which(isDE_lrt)) > 2)
{
    perform_hiearchical(counts_matrix_logCPM, isDE_lrt, cell_vec, groupnames_vec, output_folder, "lrt", "CPM_log");
    perform_hiearchical(counts_matrix_rpkm, isDE_lrt, cell_vec, groupnames_vec, output_folder, "lrt", "RPKM_log");
}

if(length(which(isDE_qlf)) > 2)
{
    perform_hiearchical(counts_matrix_logCPM, isDE_qlf, cell_vec, groupnames_vec, output_folder, "qlf", "CPM_log");
    perform_hiearchical(counts_matrix_rpkm, isDE_qlf, cell_vec, groupnames_vec, output_folder, "qlf", "RPKM_log");
}


cat("\n writing tables\n");

table_exact_up <- return_table_up_or_down(exact_nb_signif_twofolds, 1, exactest_results); 
table_exact_down <- return_table_up_or_down(exact_nb_signif_twofolds, -1, exactest_results); 
table_lrt_up <- return_table_up_or_down(results_lrt_nb_signif_twofolds, 1, results_lrt);
table_lrt_down <- return_table_up_or_down(results_lrt_nb_signif_twofolds, -1, results_lrt);
table_qlf_up <- return_table_up_or_down(results_qlf_nb_signif_twofolds, 1, results_qlf);
table_qlf_down <- return_table_up_or_down(results_qlf_nb_signif_twofolds, -1, results_qlf);



if(!is.na(table_exact_up[[1]]) && !is.na (table_exact_down[[1]]) && (length(which(isDE_exact)) != (nrow(table_exact_up) + nrow(table_exact_down))))
{
    stop("\n Problem in constructing the tables exact test\n");
}

if(!is.na(table_lrt_up[[1]]) && !is.na (table_lrt_down[[1]]) && (length(which(isDE_lrt)) != (nrow(table_lrt_up) + nrow(table_lrt_down))))
{
    stop("\n Problem in constructing the tables lrt\n");
}

if(!is.na(table_qlf_up[[1]]) && !is.na (table_qlf_down[[1]]) && (length(which(isDE_qlf)) != (nrow(table_qlf_up) + nrow(table_qlf_down))))
{
    stop("\n Problem in constructing the tables qlf\n");
}


write.table(table_exact_up, file=paste(output_folder, "table_exact_up_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);
write.table(table_exact_down, file=paste(output_folder,"table_exact_down_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);

write.table(table_lrt_up, file=paste(output_folder,"table_lrt_up_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);
write.table(table_lrt_down, file=paste(output_folder,"table_lrt_down_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);

write.table(table_qlf_up, file=paste(output_folder,"table_qlf_up_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);
write.table(table_qlf_down, file=paste(output_folder,"table_qlf_down_twofolds_signif.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T);

