#################
# This script takes a Position Weight Matrix (PWM) and output heatmaps and average profiles of the motifs.
# Descostes dec 2016
#################


library("seqPattern");
library("Biostrings");
library("NGSprofiling");
library("Rargs");



#################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


paramsDefinition[["--fastaFile"]] <- list(variableName="fasta_file", numeric=F, mandatory=T, description="A unique path to the fasta file to use.");
paramsDefinition[["--flankUp"]] <- list(variableName="flank_up", numeric=T, mandatory=T, description="A unique number indicating the interval before the centering zone. (no value = NULL)");
paramsDefinition[["--flankDown"]] <- list(variableName="flank_down", numeric=T, mandatory=T, description="A unique number indicating the interval after the centering zone. (no value = NULL)");
paramsDefinition[["--Nbin"]] <- list(variableName="n_bin", numeric=T, mandatory=T, description="A unique number indicating the bin size representation on the heatmap. (no value = NULL)");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="A unique path to the output folder.");
paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="A unique number of cpu to use.");
paramsDefinition[["--consensusIupacVec"]] <- list(variableName="consensus_iupac_vec", numeric=F, mandatory=T, description="Space separated iupac motifs to search in the fasta file provided.");
paramsDefinition[["--consensusToolsVec"]] <- list(variableName="consensus_tools_vec", numeric=F, mandatory=T, description="Space separated strings of tool names that were used to find the iupac codes.");
paramsDefinition[["--consensusMotifName_vec"]] <- list(variableName="consensus_motif_name_vec", numeric=F, mandatory=T, description="Vector of strings of the motif represented by the iupac codes.");
paramsDefinition[["--pwmVec"]] <- list(variableName="pwm_vec", numeric=F, mandatory=T, description="Space separated pathes to the pwm matrices of the motifs. Should have the same number of elements than the iupac vec.");
paramsDefinition[["--pwmNameVec"]] <- list(variableName="pwm_name_vec", numeric=F, mandatory=T, description="Space separated strings of tool names that were used to find the pwm matrices with motif name as prefix.");
paramsDefinition[["--minScorePwm"]] <- list(variableName="min_score_pwm", numeric=F, mandatory=T, description="A unique string giving the percentage of homology for searching with pwm matrices.");



################



################
# FUNCTION
################




sequence_ordering_by_position <- function(occurence_list, total_nb_sequences){
	
	#Ordering by motif position
	
	seq_index_ordered_position_list <- lapply(occurence_list, function(x){
				
				if(!is.null(x) && nrow(x) > 1)
				{
					ordered_matrix <- x[order(x$position),];
					sequence_number <- ordered_matrix[-which(duplicated(ordered_matrix$sequence)),1]; # Keeping only the left most hit for the motif
					return(sequence_number);
				}else{
					return(NULL);
				}
				
			})
	
	#Retrieving the index of sequences not having a motif match
	
	missing_index_list <- lapply(seq_index_ordered_position_list, function(x, total_nb_seq){
				
				if(!is.null(x))
				{
					return(which(is.na(match(1:total_nb_seq, x))));
				}else{
					return(NULL);
				}
				
				
			}, total_nb_sequences)
	
	
	index_null_seq <- which(unlist(lapply(seq_index_ordered_position_list,is.null)));
	missing_null_seq <- which(unlist(lapply(missing_index_list,is.null)));
	seq_index_for_test <- seq_index_ordered_position_list;
	seq_missing_for_test <- missing_index_list;
	
	if(length(index_null_seq) != 0)
	{
		seq_index_for_test <- seq_index_for_test[-index_null_seq];
	}
	
	if(length(missing_null_seq) != 0)
	{
		seq_missing_for_test <- seq_missing_for_test[-missing_null_seq];
	}
	
	
	if(length(unique(mapply(sum, unlist(lapply(seq_index_for_test,length)), unlist(lapply(seq_missing_for_test,length))))) != 1)
	{
		stop("Problem in position sorting\n");
	}
	
	
	#Merging the sorted and missing idexing
	
	seq_order_by_position_list <- mapply(function(x,y){
			
				return(list(c(x,y)));
			
			}, seq_index_ordered_position_list, missing_index_list);
	
	return(seq_order_by_position_list);
}


################



##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


checkingOutputFolder(output_folder);

if(length(fasta_file) != 1)
{
	stop("The input fasta file should be a unique path\n\n");
}

if(length(consensus_iupac_vec) != length(consensus_tools_vec))
{
	stop("consensus_iupac_vec and consensus_tools_vec should have the same length\n");
}

if(length(consensus_motif_name_vec) != length(consensus_iupac_vec))
{
	stop("One motif name is needed by iupac code\n");
}


if(length(pwm_vec) != length(pwm_name_vec))
{
	stop("One name for each pwm should be provided\n");
}

if(length(min_score_pwm) != 1)
{
	stop("Only one percentage for motif search with pwm can be given.\n");
}

cat("####### Preparing and sorting data #######\n\n");

cat("Reading the fasta file\n");
fasta_dnaStringSet <- readDNAStringSet(fasta_file);
total_nb_sequences <- length(fasta_dnaStringSet);
original_seq_order <- list(1:length(fasta_dnaStringSet));
names(original_seq_order) <- "original_order";

cat("Finding occurence of consensus motifs\n");
list_occurence_consensus_iupac <- getPatternOccurrenceList(regionsSeq = fasta_dnaStringSet, 
		patterns = consensus_iupac_vec, 
		seqOrder = original_seq_order[[1]], 
		useMulticore = TRUE, 
		nrCores = nb_cpu);
names(list_occurence_consensus_iupac) <- consensus_tools_vec;


png(filename=paste(output_folder, "nbmotifs-iupac.png",sep=""), width = 600, height = 600, bg = "transparent")
barplot(unlist(lapply(list_occurence_consensus_iupac,function(x){return(length(unique(x$sequence)))})), horiz=T, las=1, xlim=c(0,length(fasta_dnaStringSet)), xlab="nb of seq with motif");
dev.off();


cat("\t Ordering by position\n");
seq_order_by_position_list_consensus_iupac <- sequence_ordering_by_position(list_occurence_consensus_iupac, total_nb_sequences);



cat("Finding occurence of pwm\n");
list_occurence_pwm <- list();
for(i in 1:length(pwm_vec)) 
{
	cat("\t", i, "/", length(pwm_vec),"\n");
	pwm_mat <- read.table(pwm_vec[i], stringsAsFactors=FALSE);
	pwm_mat <- t(pwm_mat);
	rownames(pwm_mat) <- c("A", "C", "G", "T");
	
	list_occurence_pwm[[i]] <- motifScanHits(regionsSeq = fasta_dnaStringSet, 
		motifPWM = pwm_mat, 
		minScore = min_score_pwm);
}
names(list_occurence_pwm) <- pwm_name_vec;

png(filename=paste(output_folder, "nbmotifs-pwm.png",sep=""), width = 600, height = 600, bg = "transparent")
barplot(unlist(lapply(list_occurence_pwm,function(x){return(length(unique(x$sequence)))})), horiz=T, las=1, xlim=c(0,length(fasta_dnaStringSet)), xlab="nb of seq with motif");
dev.off();


cat("\t Ordering by position\n");
seq_order_by_position_list_pwm <- sequence_ordering_by_position(list_occurence_pwm, total_nb_sequences);
seq_order_by_position_list_consensus_iupac <- c(seq_order_by_position_list_consensus_iupac, original_seq_order);



cat("\n\n ####### Plotting heatmaps and average profiles #######\n\n");

length_sequence <- unique(width(fasta_dnaStringSet));

		
cat("Plotting motif with consensus sequences\n");


names_iupac_consensus_tool_vec <- names(seq_order_by_position_list_consensus_iupac); 
		
for(i in 1:length(seq_order_by_position_list_consensus_iupac)) 
{
	name_iupac_consensus_tool <- names_iupac_consensus_tool_vec[i];
	consensus_motif_name <- consensus_motif_name_vec[i];
	
	cat("\t", i, "/", length(seq_order_by_position_list_consensus_iupac), "\n");
	output_folder_consensus <- paste(output_folder, "consensus_iupac/", if(name_iupac_consensus_tool != "original_order") "ordered/", name_iupac_consensus_tool, "/", sep="");
	checkingOutputFolder(output_folder_consensus);
	
	if(!is.null(seq_order_by_position_list_consensus_iupac[[i]]))
	{
		plotPatternDensityMap(regionsSeq = fasta_dnaStringSet, 
				patterns = if(name_iupac_consensus_tool == "original_order"){
							
							if(length(which(unlist(lapply(seq_order_by_position_list_consensus_iupac,is.null)))) != 0)
							consensus_iupac_vec[-which(unlist(lapply(seq_order_by_position_list_consensus_iupac,is.null)))]
						else consensus_iupac_vec
							
						}  else consensus_iupac_vec[i], 
				seqOrder = seq_order_by_position_list_consensus_iupac[[i]],
				flankUp = flank_up, 
				flankDown = flank_down, 
				nBin = n_bin,  
				color = "cyan", 
				xTicks = NULL, 
				xTicksAt = NULL, 
				xLabel = "", 
				yTicks = NULL, 
				yTicksAt = NULL, 
				yLabel = "",
				addPatternLabel = FALSE,  
				addReferenceLine = TRUE, 
				plotColorLegend = TRUE,
				outFile = paste(output_folder_consensus, consensus_motif_name, if(name_iupac_consensus_tool != "original_order") paste("_",name_iupac_consensus_tool,sep=""), sep=""), 
				useMulticore = TRUE, 
				nrCores = nb_cpu)
		
	}else{ #plot an empty matrix
		
		png(filename= paste(output_folder_consensus, consensus_motif_name, "_",name_iupac_consensus_tool, "_", consensus_iupac_vec[i], ".png", sep=""), width = 600, height = 600, bg = "transparent")
		plot(NULL,main="",xlim= c(0, length_sequence), ylim= c(0,length(fasta_dnaStringSet)), xlab="", ylab="", axes=F);
		abline(v= round(length_sequence/2),lty=2, col="black",lwd=2);
		axis(side=1,labels=c(-round(length_sequence/2), -round(length_sequence/4), 0, round(length_sequence/4), round(length_sequence/2)),
				at=c(0, round(length_sequence/4), round(length_sequence/2), round(length_sequence/2) + round(length_sequence/4), length_sequence), lwd=0, lwd.ticks=1, cex.axis=1.4);
		box();
		dev.off();
	}
	

	if(name_iupac_consensus_tool == "original_order")
	{
		file_to_rename_vec <- paste(output_folder_consensus, consensus_motif_name, ".", consensus_iupac_vec, ".png", sep="");
		new_file_name_vec <- paste(output_folder_consensus, consensus_motif_name, "_", consensus_tools_vec, ".png", sep="");
		file.rename(from = file_to_rename_vec, to = new_file_name_vec);
		
	}
	
	
	## Plotting dinucleotide information
	
	cat("# Plotting dinucleotide and IUPAC ambiguity heatmaps\n");
	
	
	output_folder_dinucleotides <- paste(output_folder_consensus, "dinucleotides/", sep="");
	checkingOutputFolder(output_folder_dinucleotides);
	
	dinucleotide_patterns <- c("AA", "TA", "CG", "GC", "TT", "GG", "CC", "SS", "WW", "GA", "RR", "R", "AG"); 
	
	plotPatternDensityMap(regionsSeq = fasta_dnaStringSet, 
			patterns = dinucleotide_patterns, 
			seqOrder = if(!is.null(seq_order_by_position_list_consensus_iupac[[i]])) seq_order_by_position_list_consensus_iupac[[i]] else original_seq_order[[1]],
			flankUp = flank_up, 
			flankDown = flank_down, 
			nBin = n_bin,  
			color = "blue", 
			xTicks = NULL, 
			xTicksAt = NULL, 
			xLabel = "", 
			yTicks = NULL, 
			yTicksAt = NULL, 
			yLabel = "",
			addPatternLabel = FALSE, 
			addReferenceLine = TRUE, 
			plotColorLegend = TRUE,
			outFile = paste(output_folder_dinucleotides, "DinucleotideDensityMap", sep=""), 
			useMulticore = TRUE, 
			nrCores = nb_cpu);
	
	
	#Plotting average profiles

	cat("Plotting average profiles\n");
	
	if(name_iupac_consensus_tool != "original_order")
	{
		png(filename=paste(output_folder_consensus, consensus_motif_name, "_", name_iupac_consensus_tool, "-iupac.png",sep=""), width = 600, height = 600, bg = "transparent")
		plotPatternOccurrenceAverage(regionsSeq = fasta_dnaStringSet, 
				patterns = consensus_iupac_vec[i], 
				flankUp = flank_up, 
				flankDown = flank_down, 
				smoothingWindow = 1, 
				color = "cyan",
				xLabel = "Distance to reference point (bp)", 
				yLabel = "Relative frequency",
				plotLegend = TRUE,
				useMulticore = TRUE, 
				nrCores = nb_cpu);
		dev.off();
	}	
	
	png(filename=paste(output_folder_dinucleotides, "dinucleotide_avp.png",sep=""), width = 800, height = 800, bg = "transparent")
	plotPatternOccurrenceAverage(regionsSeq = fasta_dnaStringSet, 
			patterns = dinucleotide_patterns, 
			flankUp = flank_up, 
			flankDown = flank_down, 
			smoothingWindow = 1, 
			color = rainbow(length(dinucleotide_patterns)),
			xLabel = "Distance to reference point (bp)", 
			yLabel = "Relative frequency",
			plotLegend = TRUE,
			useMulticore = TRUE, 
			nrCores = nb_cpu)
	dev.off();
	
	
}





cat("Plotting motif with PWM\n");


for(i in 1:length(seq_order_by_position_list_pwm)) 
{
	pwm_name <- pwm_name_vec[i];
	output_folder_pwm <- paste(output_folder, "pwm_matrix/ordered/", pwm_name, "/", sep="");
	checkingOutputFolder(output_folder_pwm);
	
	pwm_mat <- read.table(pwm_vec[i], stringsAsFactors=FALSE);
	pwm_mat <- t(pwm_mat);
	rownames(pwm_mat) <- c("A", "C", "G", "T");
	
	if(!is.null(seq_order_by_position_list_pwm[[i]]))
	{
		plotMotifDensityMap(regionsSeq = fasta_dnaStringSet, 
				motifPWM = pwm_mat, 
				minScore = min_score_pwm,
				seqOrder = seq_order_by_position_list_pwm[[i]], 
				flankUp = flank_up, 
				flankDown = flank_down,
				nBin = n_bin, 
				color = "green", 
				xTicks = NULL,
				xTicksAt = NULL, 
				xLabel = "", 
				yTicks = NULL, 
				yTicksAt = NULL, 
				yLabel = "",
				plotColorLegend = TRUE, 
				outFile = paste(output_folder_pwm, pwm_name, "-pwm", sep=""));
		
		plotMotifScanScores(regionsSeq = fasta_dnaStringSet, 
				motifPWM = pwm_mat, 
				seqOrder = seq_order_by_position_list_pwm[[i]],
				flankUp = flank_up, 
				flankDown = flank_down, 
				xTicks = NULL, 
				xTicksAt = NULL,
				xLabel = "", 
				yTicks = NULL, 
				yTicksAt = NULL, 
				yLabel = "", 
				plotColorLegend = TRUE, 
				outFile = paste(output_folder_pwm, pwm_name, "-scanScores.png", sep=""));
		
		plotMotifDensityMap(regionsSeq = fasta_dnaStringSet, 
				motifPWM = pwm_mat, 
				minScore = min_score_pwm,
				seqOrder = original_seq_order[[1]], 
				flankUp = flank_up, 
				flankDown = flank_down,
				nBin = n_bin, 
				color = "green", 
				xTicks = NULL,
				xTicksAt = NULL, 
				xLabel = "", 
				yTicks = NULL, 
				yTicksAt = NULL, 
				yLabel = "",
				plotColorLegend = TRUE, 
				outFile = paste(output_folder_pwm, pwm_name, "originalOrder-pwm", sep=""))
		
		
		plotMotifScanScores(regionsSeq = fasta_dnaStringSet, 
				motifPWM = pwm_mat, 
				seqOrder = original_seq_order[[1]],
				flankUp = flank_up, 
				flankDown = flank_down, 
				xTicks = NULL, 
				xTicksAt = NULL,
				xLabel = "", 
				yTicks = NULL, 
				yTicksAt = NULL, 
				yLabel = "", 
				plotColorLegend = TRUE, 
				outFile = paste(output_folder_pwm, pwm_name, "-scanScores_originalOrder.png", sep=""));
		
		
	}else{
		
		png(filename= paste(output_folder_pwm, pwm_name, "-pwm.png", sep=""), width = 600, height = 600, bg = "transparent")
		plot(NULL,main="",xlim= c(0, length_sequence), ylim= c(0,length(fasta_dnaStringSet)), xlab="", ylab="", axes=F);
		abline(v= round(length_sequence/2),lty=2, col="black",lwd=2);
		axis(side=1,labels=c(-round(length_sequence/2), -round(length_sequence/4), 0, round(length_sequence/4), round(length_sequence/2)),
				at=c(0, round(length_sequence/4), round(length_sequence/2), round(length_sequence/2) + round(length_sequence/4), length_sequence), lwd=0, lwd.ticks=1, cex.axis=1.4);
		box();
		dev.off();
		
		png(filename= paste(output_folder_pwm, pwm_name, "-scanScores.png", sep=""), width = 600, height = 600, bg = "transparent")
		plot(NULL,main="",xlim= c(0, length_sequence), ylim= c(0,length(fasta_dnaStringSet)), xlab="", ylab="", axes=F);
		abline(v= round(length_sequence/2),lty=2, col="black",lwd=2);
		axis(side=1,labels=c(-round(length_sequence/2), -round(length_sequence/4), 0, round(length_sequence/4), round(length_sequence/2)),
				at=c(0, round(length_sequence/4), round(length_sequence/2), round(length_sequence/2) + round(length_sequence/4), length_sequence), lwd=0, lwd.ticks=1, cex.axis=1.4);
		box();
		dev.off();
		
		png(filename= paste(output_folder_pwm, pwm_name, "originalOrder-pwm.png", sep=""), width = 600, height = 600, bg = "transparent")
		plot(NULL,main="",xlim= c(0, length_sequence), ylim= c(0,length(fasta_dnaStringSet)), xlab="", ylab="", axes=F);
		abline(v= round(length_sequence/2),lty=2, col="black",lwd=2);
		axis(side=1,labels=c(-round(length_sequence/2), -round(length_sequence/4), 0, round(length_sequence/4), round(length_sequence/2)),
				at=c(0, round(length_sequence/4), round(length_sequence/2), round(length_sequence/2) + round(length_sequence/4), length_sequence), lwd=0, lwd.ticks=1, cex.axis=1.4);
		box();
		dev.off();
		
		png(filename= paste(output_folder_pwm, pwm_name, "-scanScores_originalOrder.png", sep=""), width = 600, height = 600, bg = "transparent")
		plot(NULL,main="",xlim= c(0, length_sequence), ylim= c(0,length(fasta_dnaStringSet)), xlab="", ylab="", axes=F);
		abline(v= round(length_sequence/2),lty=2, col="black",lwd=2);
		axis(side=1,labels=c(-round(length_sequence/2), -round(length_sequence/4), 0, round(length_sequence/4), round(length_sequence/2)),
				at=c(0, round(length_sequence/4), round(length_sequence/2), round(length_sequence/2) + round(length_sequence/4), length_sequence), lwd=0, lwd.ticks=1, cex.axis=1.4);
		box();
		dev.off();
		
	}
	
	
	
	png(filename=paste(output_folder_pwm, pwm_name,"-pwm-avp.png",sep=""), width = 600, height = 600, bg = "transparent")
	plotMotifOccurrenceAverage(regionsSeq = fasta_dnaStringSet, 
			motifPWM = pwm_mat, 
			minScore = min_score_pwm, 
			flankUp = flank_up, 
			flankDown = flank_down, 
			smoothingWindow = 1, 
			color = "green", 
			xLabel = "Distance to reference point (bp)", 
			yLabel = "Relative frequency", 
			plotLegend = FALSE)
	dev.off();
	
	
}










