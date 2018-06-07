###############
# This script performs gene set enrichment analysis from chip-seq peaks using the package Chip-enrich.
# Descostes december 2016
###############

library("chipenrich");
library("NGSprofiling");
library("Rargs");
library("RColorBrewer");
library("ggplot2");
library("grid");

################
# PARAMETERS
################



#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--peaksFile"]] <- list(variableName="peaks_file", numeric=F, mandatory=T, description="Peaks files can be in broadPeak, narrowPeak or in bed format.");
paramsDefinition[["--outputPrefix"]] <- list(variableName="output_prefix", numeric=F, mandatory=T, description="Prefix attached to all output file.");
paramsDefinition[["--outputFolderBase"]] <- list(variableName="output_folder_base", numeric=F, mandatory=T, description="Base name of the output folder.");
paramsDefinition[["--genomeVersion"]] <- list(variableName="genome_version", numeric=F, mandatory=T, description="ex: 'mm10', 'hg19'.");
paramsDefinition[["--genesSetsVector"]] <- list(variableName="genesSets_vector", numeric=F, mandatory=T, description="possible values are c('GOBP', 'GOCC', 'GOMF', 'reactome', 'biocarta_pathway', 'ctd', 'cytoband', 'drug_bank', 'ehmn_pathway_gene', 'gene_expression', 'hallmark', 'immunologic', 'kegg_pathway', 'mesh', 'metabolite', 'microrna', 'oncogenic', 'panther_pathway', 'pfam', 'protein_interaction_biogrid', 'transcription_factors').");
paramsDefinition[["--locus"]] <- list(variableName="locus", numeric=F, mandatory=T, description="Possible values are nearest_tss, nearest_gene, exon, intron, 1kb, 5kb, 10kb, 1kb_outside_upstream, 5kb_outside_upstream, 10kb_outside_upstream, 1kb_outside, 5kb_outside, and 10kb_outside.");
paramsDefinition[["--readLength"]] <- list(variableName="read_length", numeric=T, mandatory=T, description="Length of reads of the experiment and should not be higher than 100.");
paramsDefinition[["--numPeakThresholdVec"]] <- list(variableName="num_peak_threshold_vec", numeric=T, mandatory=T, description="Vector of integer giving the number of peaks needed for a gene to be considered.");
paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Nb of cpu used to run the analysis.");
paramsDefinition[["--FDR"]] <- list(variableName="FDR", numeric=T, mandatory=T, description="Threshold for the false discovery rate.");
paramsDefinition[["--mappabilityOption"]] <- list(variableName="mappability_option", numeric=F, mandatory=T, description="Boolean indicating if mappability should be taken into account.", postConversion=as.logical);




#peaks_file <- "/home/pouce/Documents/analysis/auts2/gene_ontology/CA2_peaks.broadPeak"; # Peaks files can be in broadPeak, narrowPeak or in bed format
#
#output_prefix <- "auts2_cortical";
#
#output_folder_base <- "/home/pouce/Documents/analysis/auts2/gene_ontology/auts2_peaks/"
#
#genome_version <- "mm10";
#
#genesSets_vector <- c("GOBP", "GOCC", "GOMF", "biocarta_pathway", "cytoband", "drug_bank", "ehmn_pathway_gene", "gene_expression", "kegg_pathway", "mesh", "metabolite", "mirbase", "panther_pathway", "pfam",                    
#		"protein_interaction_mimi", "reactome", "transcription_factors"); # possible values are c("GOBP", "GOCC", "GOMF", "biocarta_pathway", "cytoband", "drug_bank", "ehmn_pathway_gene", "gene_expression", "kegg_pathway", "mesh", "metabolite", "mirbase", "panther_pathway", "pfam", "protein_interaction_mimi", "reactome", "transcription_factors");    
#
#locus <- "1kb"; # Possible values are "10kb", "10kb_and_more_upstream", "1kb", "5kb", "exon", "intron", "nearest_gene", "nearest_tss"        
#
#read_length <- 46;  
#
#num_peak_threshold_vec <- c(1, 5, 10);
#
#nb_cpu <- 1;
#
#FDR <- 0.001



################




##############
# FUNCTIONS
###############

compute_ggplots <- function(input_table, outputFolder, top20=FALSE){
	
	ggplot_list <- list();			
	ggplot_list <- lapply(input_table, function(x, select_top){
				
				if(is.null(x) || nrow(x) == 0)
				{
					return(NULL);
				}else{
				
					title_graph <- paste(strsplit(unique(x[["Geneset.Type"]]), " ")[[1]], collapse="_");
					table_for_ggplot <- data.frame(number=nrow(x):1, values=-log10(x[["FDR"]]), names = x[["Description"]]);
					
					
					if(select_top)
					{
						table_for_ggplot <- table_for_ggplot[1:20,];
						x <- x[1:20,];
					}
					
					size_text <- min(c(1.5, round(2*(123/nrow(table_for_ggplot)), digits=1))); 
					p <- ggplot(data=table_for_ggplot, aes(x=number,y=values)) + geom_bar(stat="identity", fill= "firebrick4") + coord_flip() + geom_text(aes(label=names), hjust=c(rep("left", nrow(table_for_ggplot))), color="black", size= size_text) + theme(panel.background= element_blank()) + labs(title= title_graph) + ylab("-log10(FDR)") + xlab("");
					write.table(x, file=paste(outputFolder, title_graph, if(select_top) "-top20", ".txt", sep=""), sep="\t", quote=F, col.names=T, row.names=T);
				}
				
				return(p);}, top20)
	
	return(ggplot_list);
}



#Found at http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

#multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#	
#	# Make a list from the ... arguments and plotlist
#	plots <- c(list(...), plotlist)
#	
#	numPlots = length(plots)
#	
#	# If layout is NULL, then use 'cols' to determine layout
#	if (is.null(layout)) {
#		# Make the panel
#		# ncol: Number of columns of plots
#		# nrow: Number of rows needed, calculated from # of cols
#		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#				ncol = cols, nrow = ceiling(numPlots/cols))
#	}
#	
#	if (numPlots==1) {
#		print(plots[[1]])
#		
#	} else {
#		# Set up the page
#		grid.newpage()
#		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#		
#		# Make each plot, in the correct location
#		for (i in 1:numPlots) {
#			# Get the i,j matrix positions of the regions that contain this subplot
#			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#			
#			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#							layout.pos.col = matchidx$col))
#		}
#	}
#}

################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


#Checking param


if(length(peaks_file) != 1 || length(output_prefix) != 1 || length(locus) != 1)
{
	stop("Only one file of peaks, one prefix name and one locus can be submitted.\n\n");
}

check_genomes <- supported_genomes();

if(is.na(match(genome_version, check_genomes)))
{
	stop(paste(genome_version, " is not supported\n", sep=""));
}

check_genesets <- supported_genesets()

if(any(is.na(match(genesSets_vector, check_genesets$geneset))))
{
	stop("One or more gene sets defined is not supported. check supported genesets with supported_genesets() function\n\n");
}

check_gene_locus <- supported_locusdefs();

if(is.na(match(locus, check_gene_locus$locusdef)))
{
	stop("\n The locus defined is not supported\n\n");
}

if(read_length > 100)
{
	stop("\n Read length is too long and not supported by the package. See if this can be modified\n");
}

#Retrieve the closest supported gene length

read_length_vec <- supported_read_lengths();
read_length <- read_length_vec[which.min(abs(read_length - read_length_vec))]

# Load the table of peaks

peaks_table <- read.table(peaks_file, stringsAsFactors=F);



!!!!!!!!!!!!!!!!
!! The broadenrich method is now its own function, broadenrich(), instead of chipenrich(..., method = 'broadenrich', ...).
!!!!!!!!!!!!!!!
        

!!!!!!!!!!!!!!!!
        A new method for enrichment, polyenrich() is designed for gene set enrichment
of experiments where the presence of multiple peaks in a gene is accounted
for in the model. Use the polyenrich() function for this method.
!!!!!!!!!!!!!!!!!!!!

        
!!!!!!!!!!!!!!!!!!!!!!!!!
User interface for mappability has been streamlined. 'mappability' parameter
in broadenrich(), chipenrich(), and polyenrich() functions replaces the
three parameters previously used: 'use_mappability', 'mappa_file', and
'read_length'. The unified 'mappability' parameter can be 'NULL', a file path,
or a string indicating the read length for mappability, e.g. '24'.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
o A formerly hidden API for randomizations to assess Type I Error rates for
data sets is now exposed to the user. Each of the enrich functions has a
'randomization' parameter. See documentation and vignette for details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
Input files are read according to their file extension. Supported extensions
are bed, gff3, wig, bedGraph, narrowPeak, and broadPeak. Arbitrary extensions
are also supported, but there can be no header, and the first three columns
must be chr, start, and end.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Calling the broadenrich method with chipenrich(..., method = 'broadenrich', ...)
is no longer valid. Instead, use broadenrich().
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
for(method_enrichment in c("chipenrich", "broadenrich", "fet"))
{
	
	for(num_peak_threshold in num_peak_threshold_vec)
	{
		cat("#### Performing analysis for: ", method_enrichment, " ", mappability_option, " ", num_peak_threshold, "\n\n\n");
		output_folder <- paste(output_folder_base, method_enrichment, "/", if(mappability_option) "with_mappability/" else "no_mappability/", "min_nb_peak_", num_peak_threshold, "/", locus, "/", sep="");
		checkingOutputFolder(output_folder);
		
		result <- chipenrich(peaks_file, 
				out_name = output_prefix, 
				out_path = output_folder, 
				genome = genome_version,
				genesets = genesSets_vector, 
				locusdef = locus, 
				method = method_enrichment, 
				use_mappability = mappability_option, 
				read_length = read_length, 
				max_geneset_size = 2000, 
				num_peak_threshold = num_peak_threshold, 
				n_cores= nb_cpu);
		
		save(result, file=paste(output_folder, "results.Rdat", sep=""));
		
		
		
		## Output a barplot of the number of peaks found in genes
		
		png(filename=paste(output_folder, "intergenic_vs_intragenic.png",sep=""), width = 600, height = 600, bg = "transparent")
		barplot(c(nrow(result$peaks), nrow(peaks_table)-nrow(result$peaks)), names.arg=c("intragenic","intergenic"), ylab="Nb of peaks", col= brewer.pal(3, "Set3"), main= paste(nrow(result$peaks), " intragenic peaks", sep=""));
		dev.off();
		
		## Separate the result table per category
		
		table_per_category <- split(result$results, result$results[["Geneset.Type"]]);
		
		# Filtering each category according to the FDR and odds.Ratio 
		
		table_per_category_filtered <- list();
		
		table_per_category_filtered <-  lapply(table_per_category, function(x, FDR_value, output_folder){
					
					index_to_retrieve <- which(x[["FDR"]] <= FDR_value & x[["Odds.Ratio"]] >= 1 & x[["Status"]] == "enriched");
					
					if(length(index_to_retrieve) == 0){
						return(NULL);
					}else{
						
						category_name <- paste(strsplit(unique(x[["Geneset.Type"]]), " ")[[1]], collapse="_");
						
						png(filename=paste(output_folder, category_name,"_genesets_coverage.png",sep=""), width = 600, height = 600, bg = "transparent")
						hist((x[["N.Geneset.Peak.Genes"]]*100)/x[["N.Geneset.Genes"]], xlab= "Percentage of genes covered in gene set", main= category_name);
						dev.off();
						
						new_table <- x[index_to_retrieve,];
						
						return(new_table);
					}
					
				}, FDR, output_folder);
		
		
		# Generating the pdf with results presented as barplots and output tables
		
		ggplot_list <- compute_ggplots(table_per_category_filtered, output_folder);
		ggplot_list_top20 <- compute_ggplots(table_per_category_filtered, output_folder, TRUE);
		
		
		pdf(file= paste(output_folder, output_prefix, ".pdf", sep=""), onefile=TRUE, paper="a4", colormodel="cmyk", title="GSEA_analysis", height=20,width=10)			
		
		for(i in 1:length(ggplot_list))
		{
			if(!is.null(ggplot_list[[i]])) print(ggplot_list[[i]]);
		}
		dev.off();
		
		
		pdf(file= paste(output_folder, output_prefix, "-top20.pdf", sep=""), onefile=TRUE, paper="a4", colormodel="cmyk", title="GSEA_analysis", height=20,width=10)			
		for(i in 1:length(ggplot_list_top20))
		{
			if(!is.null(ggplot_list_top20[[i]])) print(ggplot_list_top20[[i]]);
		}
		dev.off();
		
		cat("Done\n\n");
		
	}
}

		


