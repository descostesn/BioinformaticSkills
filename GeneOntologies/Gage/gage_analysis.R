###############
# This script performs GO and Gene Set Enrichment Analysis by using the GAGE package.
# Descostes Nov 2016
###############


library("gage");
library("gageData");
library("pathview");
library("AnnotationFuncs");
library("ggplot2");
library("grid");
library("NGSprofiling");
library("Rargs");

################
# PARAMETERS
################


# Define the required options
paramsDefinition=list();

paramsDefinition[["--speciesName"]]=list(variableName="species_name", numeric=F, mandatory=T, description="Should be 'mouse' or 'human'.");
paramsDefinition[["--downregulatedFile"]]=list(variableName="downregulated_file", numeric=F, mandatory=T, description="Table output by DESeq2 of downregulated genes.");
paramsDefinition[["--upregulatedFile"]]=list(variableName="upregulated_file", numeric=F, mandatory=T, description="Table output by DESeq2 of upregulated genes.");
paramsDefinition[["--outputFileName"]]=list(variableName="output_file_name", numeric=F, mandatory=T, description="Path/filename without extension.");
paramsDefinition[["--inputType"]]=list(variableName="input_type", numeric=F, mandatory=T, description="Values should be 'DESeq2' or 'edgeR'.");
paramsDefinition[["--nbCols"]]=list(variableName="nb_cols", numeric=T, mandatory=T, description="Nb of columns in the output graph.");


#Optional arguments

paramsDefinition[["--msigDBFileVec"]]=list(variableName="msigDB_file_vec", numeric=F, mandatory=F, description="File pathes to MsigDB files in GMT format.", default=NA);
paramsDefinition[["--msigDBNameVec"]]=list(variableName="msigDB_name_vec", numeric=F, mandatory=F, description="Exp names of the MsigDB files.", default=NA);








#species_name <- "human";
#downregulated_file <- "/home/pouce/Documents/analysis/Jiaray/test_geneset/downregulated_significantPval_2fold.txt"
#upregulated_file <- "/home/pouce/Documents/analysis/Jiaray/test_geneset/upregulated_significantPval_2fold.txt"
#output_file_name <- "/home/pouce/Documents/analysis/Jiaray/test_geneset/PMGSK126_vs_MGSK126_D0_withreplicates_no_batchCorrection_significantPval_2fold"
#input_type <- "DESeq2";
#
#msigDB_file_vec <- c("/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_biocarta.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_cancerModules.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_canonicalPathwaysFromSeveralDB.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_chemicalAndGeneticPerturbations.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_oncogenicSignatures.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_reactome.gmt",
#		"/home/pouce/Documents/analysis/Jiaray/msigdb_cancer/genesets_hsapiens_TFtargets.gmt");
#
#msigDB_name_vec <- c("biocarta",
#		"cancerModules",
#		"canonicalPathwaysFromSeveralDB",
#		"chemicalAndGeneticPerturbations",
#		"oncogenicSignatures",
#		"reactome",
#		"TFtargets");
#
#nb_cols <- 2;

################


################
# FUNCTION
################


perform_GSEA <- function(foldchange_vector, kegg_genesets, go_genesets, go_biological_process, go_molecular_function, go_cellular_compartments, msigDB_list, reference, sample, use_stouffer,
		use_fold, rank_test, saa_test)
{
	result <- list();
	
	result[["kegg.gs"]] <- gage(foldchange_vector, gsets= kegg_genesets, ref= reference, samp= sample, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
	result[["kegg.2d.gs"]] <- gage(foldchange_vector, gsets= kegg_genesets, ref= reference, samp= sample, same.dir=FALSE, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
	result[["go.bp"]] <- gage(foldchange_vector, gsets= go_biological_process, ref= reference, samp= sample, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
	result[["go.mf"]] <- gage(foldchange_vector, gsets= go_molecular_function, ref= reference, samp= sample, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
	result[["go.cc"]] <- gage(foldchange_vector, gsets= go_cellular_compartments, ref= reference, samp= sample, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
	
	if(!is.na(msigDB_list[1]))
	{
		for(i in 1:length(msigDB_list))
		{
			result[[names(msigDB_list)[i]]] <- gage(foldchange_vector, gsets= msigDB_list[[i]], ref= reference, samp= sample, use.stouffer= use_stouffer, use.fold= use_fold, rank.test= rank_test, saaTest= saa_test);
		}
		
	}
	
	return(result);
};



#Found at http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}




compute_values_and_color <- function(GSEA_list, ismsig)
{
	
	if(!is.null(GSEA_list$greater))
	{
		sig.up <- GSEA_list$greater[, "q.val"] < 0.1 & !is.na(GSEA_list$greater[, "q.val"]);
		
		if(ismsig) names_up <- rownames(GSEA_list$greater)[sig.up] else names_up <- unlist(lapply(strsplit(rownames(GSEA_list$greater)[sig.up], " "), function(x){return(paste(x[-1], collapse=" "))}));
		qval_up <- GSEA_list$greater[,"q.val"][sig.up];
		index_sort_up <- order(qval_up);
		qval_up_log10 <- -log10(qval_up[index_sort_up]);
		names_up <- names_up[index_sort_up];
	}else{
		qval_up_log10 <- NULL;
		names_up <- NULL;
	}
	
	if(!is.null(GSEA_list$less))
	{
		sig.down <- GSEA_list$less[, "q.val"] < 0.1 & !is.na(GSEA_list$less[, "q.val"]);
		
		if(ismsig) names_down <- rownames(GSEA_list$less)[sig.down] else names_down <- unlist(lapply(strsplit(rownames(GSEA_list$less)[sig.down], " "), function(x){return(paste(x[-1], collapse=" "))}));
		qval_down <- GSEA_list$less[,"q.val"][sig.down];
		index_sort_down <- order(qval_down);
		qval_down_log10 <- -log10(qval_down[index_sort_down]);
		names_down <- names_down[index_sort_down];
	}else{
		qval_down_log10 <- NULL;
		names_down <- NULL;
	}
	
	if(!(length(qval_up_log10) == 0 && length(qval_down_log10) == 0))
	{
		values_to_plot <- data.frame(number= length(c(qval_up_log10, qval_down_log10*-1)):1, values=c(qval_up_log10, qval_down_log10*-1), names=c(names_up, names_down))
		color_vec <- c(rep("firebrick4", length(qval_down_log10)), rep("royalblue3", length(qval_up_log10)));
		
		return(list(values_to_plot, color_vec));
	}else{
		
		return(NULL);
	}
	
	
}


compute_ggplot <- function(GSEA_list, title_string, top10, output_folder, ismsig=FALSE)
{
	result <- compute_values_and_color(GSEA_list, ismsig);
	
	if(!is.null(nrow(result[[1]])) && (top10 && nrow(result[[1]]) > 10))
	{
		positive_values <- which(result[[1]][["values"]] > 0);
		negative_values <- which(result[[1]][["values"]] < 0);
		
		index_max_positive <- min(10, length(positive_values));
		index_max_negative <- min(10, length(negative_values));
		
		indexes_to_retrieve <- c(if(index_max_positive != 0) positive_values[1:index_max_positive], if(index_max_negative != 0) negative_values[1:index_max_negative]);
		
		if(length(indexes_to_retrieve) != 0)
		{
			result[[1]] <- result[[1]][indexes_to_retrieve, ];
			result[[1]][["names"]] <- as.character(result[[1]][["names"]]);
			result[[1]][["number"]] <- length(indexes_to_retrieve):1;
			result[[2]] <- c(rep("firebrick4", index_max_negative), rep("royalblue3", index_max_positive));
		}	
	}
	
	title_graph <- title_string;
	size_text <- min(c(1.5, round(1.5*(40/length(result[[1]][["values"]])), digits=1)));
	
	output_name <- sub(" ", "_", title_string);
	if(top10) output_name <- paste(output_name, "-top10", sep="");
	write.table(result[[1]], file=paste(output_folder, output_name, ".txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE);
	
	if(!is.null(result))
	{
		p <- ggplot(data=result[[1]], aes(x=number,y=values)) + geom_bar(stat="identity", fill= result[[2]]) + coord_flip() + geom_text(aes(label=names), hjust=c(rep("left", length(which(result[[1]]["values"] >= 0))),rep("right",length(which(result[[1]]["values"] < 0)))), color="black", size= size_text) + scale_x_continuous(breaks=NULL, name="") + scale_y_continuous(name="-log10(q-value)", limits=c(min(result[[1]]["values"])-5, max(result[[1]]["values"])+5), breaks=c(-1,0,1)) + theme(panel.background= element_blank()) + labs(title= title_graph);
	}else{
		p <- ggplot(NULL) + labs(title= title_graph);
	}
	
	return(p);
}



make_ggplot_panel <- function(GSEA_object, method_name, msigDB_list, nb_cols, output_folder, top10 = FALSE)
{
	p1 <- compute_ggplot(GSEA_object[[method_name]][["kegg.gs"]], paste(method_name, ": Kegg.gs", sep=""), top10, output_folder);
	p2 <- compute_ggplot(GSEA_object[[method_name]][["kegg.2d.gs"]], paste(method_name, ": Kegg.2d.gs", sep=""), top10, output_folder);
	p3 <- compute_ggplot(GSEA_object[[method_name]][["go.bp"]], paste(method_name, ": go.bp", sep=""), top10, output_folder);
	p4 <- compute_ggplot(GSEA_object[[method_name]][["go.mf"]], paste(method_name, ": go.mf", sep=""), top10, output_folder);
	p5 <- compute_ggplot(GSEA_object[[method_name]][["go.cc"]], paste(method_name, ": go.cc", sep=""), top10, output_folder);
	
	ggsave(paste(method_name, "kegg.gs", if(top10) "_top10", ".png", sep=""), plot = p1, device = "png", path = output_folder);
	ggsave(paste(method_name, "kegg.2d.gs.", if(top10) "_top10", ".png", sep=""), plot = p2, device = "png", path = output_folder);
	ggsave(paste(method_name, "go.bp", if(top10) "_top10", ".png", sep=""), plot = p3, device = "png", path = output_folder);
	ggsave(paste(method_name, "go.mf", if(top10) "_top10", ".png", sep=""), plot = p4, device = "png", path = output_folder);
	ggsave(paste(method_name, "go.cc", if(top10) "_top10", ".png", sep=""), plot = p5, device = "png", path = output_folder);
	
	
	if(!is.na(msigDB_list[1]))
	{
		result_p_msig <- list();
				
		for(i in 1:length(msigDB_list))
		{
			result_p_msig[[i]] <- compute_ggplot(GSEA_object[[method_name]][[names(msigDB_list)[i]]], paste(method_name, ": ", names(msigDB_list)[i], sep=""), top10, output_folder, TRUE);
			ggsave(paste(method_name, names(msigDB_list)[i], if(top10) "_top10", sep="_"), plot = result_p_msig[[i]], device = "png", path = output_folder);
		}

		plotlist= result_p_msig;
	}
	
	multiplot(p1, p2, p3, p4, p5, if(!is.na(msigDB_list[1])) plotlist= result_p_msig, cols= nb_cols);	
}

################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


checkingOutputFolder(dirname(output_file_name));


if(species_name != "human" && species_name != "mouse")
{
	stop("\n The only supported species are mouse and human\n");
}

if(input_type != "DESeq2" && input_type != "edgeR")
{
	stop("\n The input type should be DESeq2 or edgeR exactly\n");
}

cat("Loading data\n");

if(species_name == "human")
{
	library(org.Hs.eg.db);
	
	kg.hsa=kegg.gsets("human");
	kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx];
	
	go.gs=go.gsets(species="human")
	go.bp=go.gs$go.sets[go.gs$go.subs$BP]
	go.mf=go.gs$go.sets[go.gs$go.subs$MF]
	go.cc=go.gs$go.sets[go.gs$go.subs$CC]
	
}else{
	
	library(org.Mm.eg.db);
	
	kg.mm=kegg.gsets("mouse");
	kegg.gs=kg.mm$kg.sets[kg.mm$sigmet.idx];
	
	go.gs=go.gsets(species="mouse");
	go.bp=go.gs$go.sets[go.gs$go.subs$BP]
	go.mf=go.gs$go.sets[go.gs$go.subs$MF]
	go.cc=go.gs$go.sets[go.gs$go.subs$CC]
	
}


#Loading the msigDB files if defined

if(!is.na(msigDB_file_vec[1]))
{
	cat("Loading the MSigDB files\n");
	
	data(egSymb);
	
	msigDB_list <- list();
	
	for(i in 1:length(msigDB_file_vec))
	{
		current.gs <- readList(msigDB_file_vec[i]);
		current.gs.sym <- lapply(current.gs, sym2eg);
		msigDB_list[[msigDB_name_vec[i]]] <- current.gs.sym;
		
		# Removing names which value is NA
		remove_name_index <- which(is.na(names(msigDB_list[[msigDB_name_vec[i]]])));
		
		if(length(remove_name_index) != 0)
		{
			msigDB_list[[msigDB_name_vec[i]]] <- msigDB_list[[msigDB_name_vec[i]]][-remove_name_index];
		}
		
		# Removing for each element the values equal to NA
		
		msigDB_list[[msigDB_name_vec[[i]]]] <- lapply(msigDB_list[[msigDB_name_vec[[i]]]], function(x){
					
					index_na <- which(is.na(x));
					percent_remove <- (length(index_na)*100)/length(x);
					
					if(percent_remove > 50)
					{
						return(NULL); #Discarding the gene set
					}
					
					
					if(length(index_na) > 0) return(x[-index_na]) else return(x);	
				})
		
		#If gene set contains no element, remove the gene set
		
		length_vec <- unlist(lapply(msigDB_list[[msigDB_name_vec[[i]]]],length));
		index_remove <- which(length_vec == 0);
		
		if(length(index_remove) != 0)
		{
			msigDB_list[[msigDB_name_vec[i]]] <- msigDB_list[[msigDB_name_vec[i]]][-index_remove];
		}
	}
	
}



downregulated <- read.table(downregulated_file, stringsAsFactors=TRUE, header=TRUE, row.names=1);
upregulated <- read.table(upregulated_file, stringsAsFactors=TRUE, header=TRUE, row.names=1);


genes_table <- rbind(downregulated, upregulated);

if(nrow(genes_table) <= 2)
{
	stop("\n Your input files are empty or contain less than two lines, no analysis performed\n\n");
}

#Converting to entrez ID

entrezID_vec <- translate(rownames(genes_table), if(species_name == "human") org.Hs.egREFSEQ2EG else org.Mm.egREFSEQ2EG, remove.missing= FALSE);

if(length(which(entrezID_vec == "NA")) != 0)
{
	entrezID_vec <- entrezID_vec[-which(entrezID_vec == "NA")];
}

if(is.unsorted(match(names(entrezID_vec), rownames(genes_table))))
{
	cat("Reordering entrezID\n");
	entrezID_vec <- entrezID_vec[match(names(entrezID_vec), rownames(genes_table))];
}

index_to_filter <- match(names(entrezID_vec), rownames(genes_table));
genes_table <- genes_table[index_to_filter,];

		
		
#Removing duplicated annotations (not optimal but no other choices when converting to entrezID)

indexes_to_remove <- which(duplicated(entrezID_vec));

if(length(indexes_to_remove) != 0)
{
	cat("\t The number of duplicated annotations due to the conversion to entrezID is: ", length(indexes_to_remove),"/", nrow(genes_table), "\n");
	genes_table <- genes_table[-indexes_to_remove,];
	entrezID_vec <- entrezID_vec[-indexes_to_remove];
}


rownames(genes_table) <- entrezID_vec;
exp.fc <- if(input_type == "DESeq2") genes_table[,2] else genes_table[,1];
names(exp.fc) <- rownames(genes_table);


cat("Performing the gene set enrichment analysis\n");


GSEA_result <- list();

cat("\t Using GAGE t-test only\n");
GSEA_result[["with_Gage_tTest"]] <- perform_GSEA(exp.fc, kegg.gs, go.gs, go.bp, go.mf, go.cc, msigDB_list, NULL, NULL, TRUE, FALSE, FALSE, gs.tTest);

cat("\t Using GAGE without stouffer\n");
GSEA_result[["with_Gage_noStouffer"]] <- perform_GSEA(exp.fc, kegg.gs, go.gs, go.bp, go.mf, go.cc, msigDB_list, NULL, NULL, FALSE, TRUE, FALSE, gs.tTest);

cat("\t Using  Mann Whitney U tests\n");
GSEA_result[["with_Gage_MWU"]] <- perform_GSEA(exp.fc, kegg.gs, go.gs, go.bp, go.mf, go.cc, msigDB_list, NULL, NULL, TRUE, TRUE, TRUE, gs.tTest);

cat("\t Using the non-parametric Kolmogorov-Smirnov tests\n")
GSEA_result[["with_Gage_KST"]] <- perform_GSEA(exp.fc, kegg.gs, go.gs, go.bp, go.mf, go.cc, msigDB_list, NULL, NULL, TRUE, TRUE, FALSE, gs.KSTest);

cat("\t Using PAGE\n");
GSEA_result[["with_PAGE"]] <- perform_GSEA(exp.fc, kegg.gs, go.gs, go.bp, go.mf, go.cc, msigDB_list, NULL, NULL, TRUE, TRUE, FALSE, gs.zTest);



cat("Formatting results\n");



pdf(file= paste(output_file_name, ".pdf", sep=""), onefile=TRUE, paper="a4", colormodel="cmyk", title="GSEA_analysis", height=20,width=10)

make_ggplot_panel(GSEA_result, "with_Gage_noStouffer", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""));
make_ggplot_panel(GSEA_result, "with_Gage_tTest", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""));
make_ggplot_panel(GSEA_result, "with_Gage_MWU", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""));
make_ggplot_panel(GSEA_result, "with_Gage_KST", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""));
make_ggplot_panel(GSEA_result, "with_PAGE", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""));

dev.off();


pdf(file= paste(output_file_name, "top10.pdf", sep=""), onefile=TRUE, paper="a4", colormodel="cmyk", title="GSEA_analysis", height=20,width=10)

make_ggplot_panel(GSEA_result, "with_Gage_noStouffer", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""), TRUE);
make_ggplot_panel(GSEA_result, "with_Gage_tTest", msigDB_list, nb_cols, paste(dirname(output_file_name), "/", sep=""), TRUE);
make_ggplot_panel(GSEA_result, "with_Gage_MWU", msigDB_list, nb_cols,  paste(dirname(output_file_name), "/", sep=""), TRUE);
make_ggplot_panel(GSEA_result, "with_Gage_KST", msigDB_list, nb_cols,  paste(dirname(output_file_name), "/", sep=""), TRUE);
make_ggplot_panel(GSEA_result, "with_PAGE", msigDB_list, nb_cols,  paste(dirname(output_file_name), "/", sep=""), TRUE);

dev.off();





