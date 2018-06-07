##############
# This script performs gene ontologies from gene lists with cluster profiler. It enables also the comparison of ontologies between different samples.
# Descostes march 2017
##############



library("clusterProfiler");
library("ReactomePA");
library("RDAVIDWebService");
library("qusage");
library("NGSprofiling");
library("topGO");
library("Rgraphviz");
library("Rargs");


################
# PARAMETERS
################

# Define the required options
paramsDefinition=list();

paramsDefinition[["--gffVec"]]=list(variableName="gff_vec", numeric=F, mandatory=T, description="Vector of gff files containing the annotations of the different categories.");
paramsDefinition[["--expnamesVec"]]=list(variableName="expnames_vec", numeric=F, mandatory=T, description="Vector of expnames corresponding to the gff files.");
paramsDefinition[["--speciesName"]]=list(variableName="species_name", numeric=F, mandatory=T, description="Single string giving the species name. Should be mouse or human.");
paramsDefinition[["--outputFolder"]]=list(variableName="output_folder", numeric=F, mandatory=T, description="Single string giving the path to the output folder.");
paramsDefinition[["--outputFormat"]]=list(variableName="output_format", numeric=F, mandatory=T, description="Should be png, pdf or ps.");

#optional
paramsDefinition[["--gmtFileVec"]]=list(variableName="gmt_file_vec", numeric=F, mandatory=F, description="Vector of gmt files (not yet supported).", default=NA);






#gff_vec <- c("/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/genes_rnaseq_june2017/no_batchcorrections/WTDay0ToDay5up_dKOcardioD5upAndDown/peaks_per_circle/dKO-down_original.gff",
#        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/genes_rnaseq_june2017/no_batchcorrections/WTDay0ToDay5up_dKOcardioD5upAndDown/peaks_per_circle/WT-Day0Day5-up_dKO-down_original.gff", 
#        "/ifs/home/descon01/analysis/fact_ledgf/vennDiagrams/genes_rnaseq_june2017/no_batchcorrections/WTDay0ToDay5up_dKOcardioD5upAndDown/peaks_per_circle/WT-Day0Day5-up_original.gff");
#
#expnames_vec <- c("dKO-down", "common", "Day0-5-up")
#
#species_name <- "mouse";
#output_folder <- c("/ifs/home/descon01/analysis/fact_ledgf/gene_ontology/cluster_profiler/venndiagram/WTDay0ToDay5up_dKOcardioD5upAndDown/");
#output_format <- "png"
#gmt_file_vec <- NA

################




################
# FUNCTION
################


perform_barplot_and_write_table <- function(result, outputFolder, methodName, output_format){
	
	#plotting the top10
	top10 <- barplot(result, drop=TRUE, showCategory=10);
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-barplot_top10.png",sep=""), plot = top10, device = "png", path = outputFolder);
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-barplot_top10.ps",sep=""), plot = top10, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-barplot_top10.pdf",sep=""), plot = top10, device = "pdf", path = outputFolder);
 }
	
	
	#plotting the top 20
	top20 <- barplot(result, drop=TRUE, showCategory=20);
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-barplot_top20.png",sep=""), plot = top20, device = "png", path = outputFolder);
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-barplot_top20.ps",sep=""), plot = top20, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-barplot_top20.pdf",sep=""), plot = top20, device = "pdf", path = outputFolder);
 }
	
	#plotting the top 50
	top50 <- barplot(result, drop=TRUE, showCategory=50);
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-barplot_top50.png",sep=""), plot = top50, device = "png", path = outputFolder);
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-barplot_top50.ps",sep=""), plot = top50, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-barplot_top50.pdf",sep=""), plot = top50, device = "pdf", path = outputFolder);
 }
	
	#plotting all terms
	all <- barplot(result, drop=TRUE, showCategory=nrow(result));
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-barplot_all.png",sep=""), plot = all, device = "png", path = outputFolder);
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-barplot_all.ps",sep=""), plot = all, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-barplot_all.pdf",sep=""), plot = all, device = "pdf", path = outputFolder);
 }
		
	cat("\t\t\t Writing table of results\n");
	write.table(as.data.frame(result), file=paste(outputFolder, methodName, "-table.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE);
}


perform_dotplot <- function(result, outputFolder, methodName, output_format)
{
	top10 <- dotplot(result, x = "geneRatio", colorBy = "p.adjust", showCategory = 10, font.size = 12, title = "")
	top20 <- dotplot(result, x = "geneRatio", colorBy = "p.adjust", showCategory = 20, font.size = 10, title = "")
	top50 <- dotplot(result, x = "geneRatio", colorBy = "p.adjust", showCategory = 50, font.size = 8, title = "")
	all <- dotplot(result, x = "geneRatio", colorBy = "p.adjust", showCategory = nrow(result), font.size = 5, title = "")
	
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-dotplot_top10.png",sep=""), plot = top10, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top20.png",sep=""), plot = top20, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top50.png",sep=""), plot = top50, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_all.png",sep=""), plot = all, device = "png", path = outputFolder);
		
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-dotplot_top10.ps",sep=""), plot = top10, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top20.ps",sep=""), plot = top20, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top50.ps",sep=""), plot = top50, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_all.ps",sep=""), plot = all, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-dotplot_top10.pdf",sep=""), plot = top10, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_top20.pdf",sep=""), plot = top20, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_top50.pdf",sep=""), plot = top50, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_all.pdf",sep=""), plot = all, device = "pdf", path = outputFolder);
 }
}



perform_enrichMap <- function(result, outputFolder, methodName, output_format)
{
	if(nrow(result) != 0)
	{
		if(output_format == "png")
		{
			png(filename=paste(outputFolder, methodName,"-enrichHeatmap.png",sep=""), width = 1000, height = 1000, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(outputFolder, methodName,"-enrichHeatmap.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(outputFolder, methodName,"-enrichHeatmap.pdf",sep=""), width=10, height=10)
  }
		enrich_map <- enrichMap(result, vertex.label.font = 3);
		dev.off();
	}
}




perform_cnetplot <- function(result, outputFolder, methodName, output_format){
	
	if(nrow(result) != 0)
	{
		if(output_format == "png")
		{
			png(filename=paste(outputFolder, methodName, "-cnetplot.png",sep=""), width = 1000, height = 1000, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(outputFolder, methodName, "-cnetplot.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(outputFolder, methodName, "-cnetplot.pdf",sep=""), width=10, height=10)
  }
		cnetplot(result, showCategory = 5, categorySize="pvalue")
		dev.off();
	}
}


perform_GOgraph <- function(result, outputFolder, methodName, output_format)
{
	if(nrow(result) != 0)
	{
		if(output_format == "png")
		{
			png(filename=paste(outputFolder, methodName, "-GOgraph.png",sep=""), width = 1000, height = 1000, bg = "transparent")
		}else if(output_format == "ps"){
			cairo_ps(filename=paste(outputFolder, methodName, "-GOgraph.ps",sep=""), width = 7, height = 7, bg = "transparent");
		}else{
      pdf(file=paste(outputFolder, methodName, "-GOgraph.pdf",sep=""), width=10, height=10)
  }
		plotGOgraph(result);
		dev.off();
	}
}


perform_dotplot_comparison <- function(result, outputFolder, methodName, output_format)
{
	top10 <- dotplot(result,colorBy = "p.adjust", showCategory = 10, font.size = 10, title = "")
	top20 <- dotplot(result, colorBy = "p.adjust", showCategory = 20, font.size = 8, title = "")
	top50 <- dotplot(result, colorBy = "p.adjust", showCategory = 50, font.size = 6, title = "")
	all <- dotplot(result, colorBy = "p.adjust", showCategory = nrow(result), font.size = 2, title = "")
	
	if(output_format == "png")
	{
		ggsave(paste(methodName,"-dotplot_top10.png",sep=""), plot = top10, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top20.png",sep=""), plot = top20, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top50.png",sep=""), plot = top50, device = "png", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_all.png",sep=""), plot = all, device = "png", path = outputFolder);
		
	}else if(output_format == "ps"){
		ggsave(paste(methodName,"-dotplot_top10.ps",sep=""), plot = top10, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top20.ps",sep=""), plot = top20, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_top50.ps",sep=""), plot = top50, device = "ps", path = outputFolder);
		ggsave(paste(methodName,"-dotplot_all.ps",sep=""), plot = all, device = "ps", path = outputFolder);
	}else{
     ggsave(paste(methodName,"-dotplot_top10.pdf",sep=""), plot = top10, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_top20.pdf",sep=""), plot = top20, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_top50.pdf",sep=""), plot = top50, device = "pdf", path = outputFolder);
     ggsave(paste(methodName,"-dotplot_all.pdf",sep=""), plot = all, device = "pdf", path = outputFolder);
 }
		
	cat("\t\t\t Writing table of results\n");
	write.table(as.data.frame(result), file=paste(outputFolder, methodName, "-table.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE);
	
}

################



##############
# MAIN
##############

# Retreives the parameters
getParams(paramsDefinition);


if(species_name != "human" && species_name != "mouse")
{
	stop("\n The only supported species are mouse and human\n");
}

if(length(gmt_file_vec) == 0)
{
	stop("\n GMT files should be provided\n");
}

if(length(gff_vec) != length(expnames_vec))
{
	stop("\n One name has to be given per gff file\n");
}

if(length(output_folder) != 1)
{
	stop("\n outputfolder should be unique\n");
}

checkingOutputFolder(output_folder);

if(species_name == "human"){

	database_name <- "org.Hs.eg.db";
	kegg_name <- "hsa";
	
}else{

	database_name <-"org.Mm.eg.db";
	kegg_name <- "mmu";
}

cat("Reading gff files and performing symbol conversion\n");

different_id_list <- list();

for(i in 1:length(gff_vec)) 
{
	cat(i, "/", length(gff_vec),"\n");
	fi <- read.table(gff_vec[i], stringsAsFactors=F);
	
	different_id_list[[i]] <- bitr(fi$V3, fromType="REFSEQ", toType=c("SYMBOL", "ENTREZID", "ENSEMBL", "UNIPROT", "ENSEMBLPROT"), OrgDb= database_name);
}


background_id_vec <- unique(unlist(lapply(different_id_list, function(x){return(unique(x$ENTREZID))})));

cat("Performing per sample investigations\n");

for(i in 1:length(different_id_list)) 
{
	cat("\t Processing ", i, "/", length(different_id_list), "\n");
	
	cat("\t\t Grouping gene ontologies\n");
	ggo_CC <- groupGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "CC", level = 2, readable = TRUE);
	ggo_BP <- groupGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "BP", level = 2, readable = TRUE);
	ggo_MF <- groupGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "MF", level = 2, readable = TRUE);
	
	cat("\t\t Computing gene ontologies enrichment\n");
	ego_CC <- enrichGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	ego_BP <- enrichGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	ego_MF <- enrichGO(gene = unique(different_id_list[[i]]$ENTREZID), OrgDb= database_name, keytype = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	
	cat("\t\t Computing KEGG enrichment\n");
	kk <- enrichKEGG(gene = unique(different_id_list[[i]]$ENTREZID), organism = kegg_name, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2, use_internal_data = FALSE);
	
	cat("\t\t Computing KEGG modules enrichment\n");
	mkk <- enrichMKEGG(gene = unique(different_id_list[[i]]$ENTREZID), organism = kegg_name, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2);
	
	cat("\t\t Computing disease related enrichment with DOSE\n");
	edo <- enrichDO(gene = unique(different_id_list[[i]]$ENTREZID), ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2, readable = TRUE);
		
	cat("\t\t Computing reactome enrichment\n");
	react <- enrichPathway(gene = unique(different_id_list[[i]]$ENTREZID), organism = species_name, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	
	
	
#	cat("\t\t Computing DAVID analysis\n");
#	if(length(unique(different_id_list[[i]]$ENTREZID)) > 3000)
#	{
#		david_first3000 <-  enrichDAVID(gene = unique(different_id_list[[i]]$ENTREZID)[1:2999], idType = "ENTREZ_GENE_ID", listType = "Gene", minGSSize = 5, maxGSSize = 1000, annotation = "GOTERM_BP_ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, 
#				species = species_name, david.user = "descon01@nyumc.org");
#	}else{
#		
#		david_first3000 <-  enrichDAVID(gene = unique(different_id_list[[i]]$ENTREZID), idType = "ENTREZ_GENE_ID", listType = "Gene", minGSSize = 5, maxGSSize = 1000, annotation = "GOTERM_BP_ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, 
#				species = species_name, david.user = "descon01@nyumc.org");
#	}
#	
	
#	cat("Perform enrichment on MSigDB gmt files\n");
#	for(j in 1:length(gmt_file_vec)) 
#	{
#		current_gmt <- read.gmt(gmt_file_vec[j]);
#		egmt <- enricher(gene = unique(different_id_list[[i]]$ENTREZID), pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2, TERM2GENE = current_gmt);
#	}
	
	
	output_folder_category <- paste(output_folder, expnames_vec[i], "/", sep="");
	checkingOutputFolder(output_folder_category);
	
	cat("\t\t Plotting barplots\n");
	perform_barplot_and_write_table(ggo_CC, output_folder_category, "ggo_CC", output_format);
	perform_barplot_and_write_table(ggo_BP, output_folder_category, "ggo_BP", output_format);
	perform_barplot_and_write_table(ggo_MF, output_folder_category, "ggo_MF", output_format);
	perform_barplot_and_write_table(ego_CC, output_folder_category, "ego_CC", output_format);
	perform_barplot_and_write_table(ego_BP, output_folder_category, "ego_BP", output_format);
	perform_barplot_and_write_table(ego_MF, output_folder_category, "ego_MF", output_format);
	perform_barplot_and_write_table(kk, output_folder_category, "kk", output_format);
	perform_barplot_and_write_table(mkk, output_folder_category, "mkk", output_format);
	perform_barplot_and_write_table(edo, output_folder_category, "edo", output_format);
	perform_barplot_and_write_table(react, output_folder_category, "react", output_format);
	#perform_barplot_and_write_table(david_first3000, output_folder_category, "davidfirst3000");
	
	cat("\t\t Plotting dotplots\n");
	perform_dotplot(ego_CC, output_folder_category, "ego_CC", output_format);
	perform_dotplot(ego_BP, output_folder_category, "ego_BP", output_format);
	perform_dotplot(ego_MF, output_folder_category, "ego_MF", output_format);
	perform_dotplot(kk, output_folder_category, "kk", output_format);
	perform_dotplot(mkk, output_folder_category, "mkk", output_format);
	perform_dotplot(edo, output_folder_category, "edo", output_format);
	perform_dotplot(react, output_folder_category, "react", output_format);
	#perform_dotplot(david_first3000, output_folder_category, "davidfirst3000");
	
	cat("\t\t Plotting maps\n");
	perform_enrichMap(ego_CC, output_folder_category, "ego_CC", output_format);
	perform_enrichMap(ego_BP, output_folder_category, "ego_BP", output_format);
	perform_enrichMap(ego_MF, output_folder_category, "ego_MF", output_format);
	perform_enrichMap(kk, output_folder_category, "kk", output_format);
	perform_enrichMap(mkk, output_folder_category, "mkk", output_format);
	perform_enrichMap(edo, output_folder_category, "edo", output_format);
	perform_enrichMap(react, output_folder_category, "react", output_format);
	#perform_enrichMap(david_first3000, output_folder_category, "davidfirst3000");
		
	cat("\t\t Plotting cNet\n");
	perform_cnetplot(ego_CC, output_folder_category, "ego_CC", output_format);
	perform_cnetplot(ego_BP, output_folder_category, "ego_BP", output_format);
	perform_cnetplot(ego_MF, output_folder_category, "ego_MF", output_format);
	perform_cnetplot(kk, output_folder_category, "kk", output_format);
	perform_cnetplot(mkk, output_folder_category, "mkk", output_format);
	perform_cnetplot(edo, output_folder_category, "edo", output_format);
	perform_cnetplot(react, output_folder_category, "react", output_format);
	#perform_cnetplot(david_first3000, output_folder_category, "davidfirst3000");
	
#	cat("\t\t Plotting cNet\n");
#	perform_GOgraph(ego_CC, output_folder_category, "ego_CC");
#	perform_GOgraph(ego_BP, output_folder_category, "ego_BP");
#	perform_GOgraph(ego_MF, output_folder_category, "ego_MF");
}
		



if(length(gff_vec) > 1)
{
	cat("Creating the input list of entrezID\n");
	different_id_list_forComparison <- lapply(different_id_list, function(x){return(unique(x$ENTREZID))});
	names(different_id_list_forComparison) <- expnames_vec;
	
	cat("Performing clusters comparison\n");
	
	cat("\t groupGO (3)\n");
	comparison_groupGO_CC <- compareCluster(geneClusters = different_id_list_forComparison, fun = "groupGO", OrgDb= database_name, keytype = "ENTREZID", ont = "CC", level = 2, readable = TRUE);
	comparison_groupGO_BP <- compareCluster(geneClusters = different_id_list_forComparison, fun = "groupGO", OrgDb= database_name, keytype = "ENTREZID", ont = "BP", level = 2, readable = TRUE);
	comparison_groupGO_MF <- compareCluster(geneClusters = different_id_list_forComparison, fun = "groupGO", OrgDb= database_name, keytype = "ENTREZID", ont = "MF", level = 2, readable = TRUE);
	
	cat("\t enrichGO (3)\n");
	comparison_enrichGO_CC <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichGO", OrgDb= database_name, keytype = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	comparison_enrichGO_BP <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichGO", OrgDb= database_name, keytype = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	comparison_enrichGO_MF <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichGO", OrgDb= database_name, keytype = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	
	cat("\t enrichKEGG\n");
	comparison_enrichKEGG <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichKEGG", organism = kegg_name, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2, use_internal_data = FALSE);
	
	cat("\t enrichDO\n");
	comparison_enrichDO <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichDO", ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, qvalueCutoff = 0.2, readable = TRUE);
	
	cat("\t enrichPathway\n");
	comparison_enrichPathway <- compareCluster(geneClusters = different_id_list_forComparison, fun = "enrichPathway", organism = species_name, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, universe = background_id_vec, minGSSize = 5, maxGSSize = 1000, readable = TRUE);
	
	
	cat("\t Output dotplots of the comparisons\n");
	
	output_folder_compare_cluster <- paste(output_folder, "compare_cluster/", sep="");
	checkingOutputFolder(output_folder_compare_cluster);
	
	perform_dotplot_comparison(comparison_groupGO_CC, output_folder_compare_cluster, "comparison_groupGO_CC", output_format);
	perform_dotplot_comparison(comparison_groupGO_BP, output_folder_compare_cluster, "comparison_groupGO_BP", output_format);
	perform_dotplot_comparison(comparison_groupGO_MF, output_folder_compare_cluster, "comparison_groupGO_MF", output_format);
	perform_dotplot_comparison(comparison_enrichGO_CC, output_folder_compare_cluster, "comparison_enrichGO_CC", output_format);
	perform_dotplot_comparison(comparison_enrichGO_BP, output_folder_compare_cluster, "comparison_enrichGO_BP", output_format);
	perform_dotplot_comparison(comparison_enrichGO_MF, output_folder_compare_cluster, "comparison_enrichGO_MF", output_format);
	perform_dotplot_comparison(comparison_enrichKEGG, output_folder_compare_cluster, "comparison_enrichKEGG", output_format);
	perform_dotplot_comparison(comparison_enrichDO, output_folder_compare_cluster, "comparison_enrichDO", output_format);
	perform_dotplot_comparison(comparison_enrichPathway, output_folder_compare_cluster, "comparison_enrichPathway", output_format);
}





