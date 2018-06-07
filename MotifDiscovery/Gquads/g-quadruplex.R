############
# This script aims at detecting G-quads in given sequences provided in fasta format.
# Descostes, ap 2016
############

library(Biostrings);
library(pqsfinder)
library(NGSprofiling);
library(parallel);
library(Rargs);


################
# PARAMETERS
################

paramsDefinition <- list();


# Required arguments
paramsDefinition[["--fastaFileVec"]] <- list(variableName="fasta_file_vec", numeric=F, mandatory=T, description="Vector of fasta files");
paramsDefinition[["--seqLength"]] <- list(variableName="seq_length", numeric=T, mandatory=T, description="Single numeric giving the number of columns of the matrix");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Folder where results will be written");
paramsDefinition[["--expname"]] <- list(variableName="expname", numeric=F, mandatory=T, description="Name of the exp used as prefix in output files");
paramsDefinition[["--nbCpu"]] <- list(variableName="nb_cpu", numeric=T, mandatory=T, description="Number of cpu");
paramsDefinition[["--objectVec"]] <- list(variableName="object_vec", numeric=F, mandatory=T, description="Objects corresponding to the fasta file");
paramsDefinition[["--sortedCoordGffVec"]] <- list(variableName="sorted_coord_gff_vec", numeric=F, mandatory=T, description="GFF of sorted genes obtained from the cluster in S");
paramsDefinition[["--spaceSizeObject"]] <- list(variableName="space_size_object", numeric=T, mandatory=T, description="Space size of the objects provided");
paramsDefinition[["--intervalUpstream"]] <- list(variableName="intervalUpstream", numeric=T, mandatory=T, description="Interval upstream TSS of the object for computing mean values");
paramsDefinition[["--intervalDownstream"]] <- list(variableName="intervalDownstream", numeric=T, mandatory=T, description="Interval downstream TSS of the object for computing mean values");
paramsDefinition[["--binsize"]] <- list(variableName="binsize", numeric=T, mandatory=T, description="Binsize of the object used");
paramsDefinition[["--profileLength"]] <- list(variableName="profile_length", numeric=T, mandatory=T, description="profile_length of the object used");


################


################
# FUNCTION
################



isG4 <- function(subject, score, start, width, loop_1,
		run_2, loop_2, run_3, loop_3, run_4) {
	r1 <- loop_1 - start
	r2 <- loop_2 - run_2
	r3 <- loop_3 - run_3
	r4 <- start + width - run_4
	
	if (!(r1 == r2 && r1 == r3 && r1 == r4))
		return(0)
	
	run_1_s <- subject[start:start+r1-1]
	run_2_s <- subject[run_2:run_2+r2-1]
	run_3_s <- subject[run_3:run_3+r3-1]
	run_4_s <- subject[run_4:run_4+r4-1]
	
	if (length(grep("^G+$", run_1_s)) && length(grep("^C+$", run_2_s)) &&
			length(grep("^G+$", run_3_s)) && length(grep("^C+$", run_4_s)))
		return(r1 * 20)
	else
		return(0)
}



exploring_gquads <- function(coord_gquad_list, output_folder, expname, mean_vec, seq_length)
{
	#Saving the list of gquads
	cat("\t Saving g-quads list\n");
	save(coord_gquad_list, file=paste(output_folder, expname, "-listgquads", sep=""));
	
	
	#Generate histogram of width
	cat("\t Generate histogram of width\n");
	width_vec <- as.numeric(unlist(lapply(coord_gquad_list, function(x){return(x[,3])})));
	cairo_ps(filename=paste(output_folder, "lengthMotif-", expname,".ps",sep=""), width = 7, height = 7, bg = "transparent")
	hist(width_vec, xlab="motif length", main=paste(expname, ": min=", round(min(width_vec)), " max=", round(max(width_vec))))
	dev.off();
	
	
	#Generate correlation between level and width
	cat("\t Generate correlation between level and width\n");
	nb_motif_vec <- unlist(lapply(coord_gquad_list,nrow)); 
	no_motif_index <- which(nb_motif_vec == 0);
	
	if(length(no_motif_index) != 0)
	{
		nb_motif_vec <- nb_motif_vec[-no_motif_index];
		mean_vec <- mean_vec[-no_motif_index];
		coord_gquad_list <- coord_gquad_list[-no_motif_index];
	}
	
	width_vec <- as.numeric(unlist(lapply(coord_gquad_list, function(x){return(x[,3])})));
	width_mean_values <- rep(mean_vec, nb_motif_vec);
	regression <- lm(width_vec~width_mean_values);
	
	cairo_ps(filename=paste(output_folder, "widthvsvalue-", expname,".ps",sep=""), width = 7, height = 7, bg = "transparent")
	plot(width_vec, width_mean_values, type="p", pch=".", xlab="motif length", ylab="binding value", main=expname);
	if(!is.na(regression$coef[1]) && !is.na(regression$coef[2]))
	{
		abline(a=regression$coef[1],b=regression$coef[2],col="red");
	}
	dev.off();
	
	
	
	#Generate histogram of scores
	cat("\t Generate histogram of scores\n");
	score_vec <- as.numeric(unlist(lapply(coord_gquad_list, function(x){return(x[,4])})));
	
	cairo_ps(filename=paste(output_folder, "score-", expname,".ps",sep=""), width = 7, height = 7, bg = "transparent")
	hist(score_vec, xlab="score", main=paste(expname, ": min=", round(min(score_vec)), " max=", round(max(score_vec))))
	dev.off();
	
	
	
	#Correlation between scores and binding value
	cat("\t Correlation between scores and binding value\n");
	regression <- lm(score_vec~width_mean_values);
	
	cairo_ps(filename=paste(output_folder, "scorevsvalue-", expname,".ps",sep=""), width = 7, height = 7, bg = "transparent")
	plot(score_vec, width_mean_values, type="p", pch=".", xlab="score", ylab="binding value", main=expname);
	if(!is.na(regression$coef[1]) && !is.na(regression$coef[2]))
	{
		abline(a=regression$coef[1],b=regression$coef[2],col="red");
	}
	dev.off();
	
	
	
	#Generating indexes to set value '1' for each sequence
	
	cat("\t Building matrix of gquads sequences\n");
	
	seq_list <- lapply(coord_gquad_list, function(x, seq_length){
				
				result <- rep(0, seq_length);
				
				if(nrow(x) != 0)
				{
					current_start <- x[,1];
					current_end <- x[,2];
					
					if(max(current_end) > seq_length) #This should not happen
					{
						current_end[which(current_end > seq_length)] <- 19999;
					}
					
					index_1 <- unlist(mapply(function(x,y){return(seq(x,y))}, current_start, current_end));
					
					result[index_1] <- 1;
				}
				
				return(result);
				
			}, seq_length);
	
	result_matrix <- do.call(rbind, seq_list);
	
	cat("Writting the files for visualization in treeview...");
	interpolatedCoordinates <- 1:ncol(result_matrix);
	annoNames <- paste("seq_", 1:nrow(result_matrix));
	writeTreeviewCluster(paste(output_folder, "matrix-", expname, sep=""), result_matrix, interpolatedCoordinates, annoNames, compress=FALSE);
	
	
	
	#Generating indexes to set score for each sequence
	
	cat("Building matrix of scores\n");
	
	seq_list_score <- lapply(coord_gquad_list, function(x, seq_length){
				
				result <- rep(0, seq_length);
				
				if(nrow(x) != 0)
				{
					current_start <- x[,1];
					current_end <- x[,2];
					current_score <- x[,4]
					
					if(max(current_end) > seq_length) #This should not happen
					{
						current_end[which(current_end > seq_length)] <- 19999;
					}
					
					index_1 <- unlist(mapply(function(x,y){return(seq(x,y))}, current_start, current_end));
					score_1 <- unlist(mapply(function(x,y,z){return(rep(z, (y-x)+1))}, current_start, current_end, current_score));
					
					result[index_1] <- score_1;
				}
				
				return(result);
				
			}, seq_length);
	
	result_matrix_score <- do.call(rbind, seq_list_score);
	
	cat("Writting the files for visualization in treeview...");
	interpolatedCoordinates_score <- 1:ncol(result_matrix_score);
	annoNames_score <- paste("seq_", 1:nrow(result_matrix_score));
	writeTreeviewCluster(paste(output_folder, "matrix-scores-", expname, sep=""), result_matrix_score, interpolatedCoordinates_score, annoNames_score, compress=FALSE);
	
	
	#Generating indexes to set width for each sequence
	
	cat("Building matrix of width\n");
	
	seq_list_width <- lapply(coord_gquad_list, function(x, seq_length){
				
				result <- rep(0, seq_length);
				
				if(nrow(x) != 0)
				{
					current_start <- x[,1];
					current_end <- x[,2];
					current_width <- x[,3]
					
					if(max(current_end) > seq_length) #This should not happen
					{
						current_end[which(current_end > seq_length)] <- 19999;
					}
					
					index_1 <- unlist(mapply(function(x,y){return(seq(x,y))}, current_start, current_end));
					width_1 <- unlist(mapply(function(x,y,z){return(rep(z, (y-x)+1))}, current_start, current_end, current_width));
					
					if(length(index_1) != length(width_1)) #should not happen
					{
						stop("\n The number of width is not equal to the number of coord\n");
					}
					
					result[index_1] <- width_1;
				}
				
				return(result);
				
			}, seq_length);
	
	result_matrix_width <- do.call(rbind, seq_list_width);
	
	cat("Writting the files for visualization in treeview...");
	interpolatedCoordinates_width <- 1:ncol(result_matrix_width);
	annoNames_width <- paste("seq_", 1:nrow(result_matrix_width));
	writeTreeviewCluster(paste(output_folder, "matrix-width-", expname, sep=""), result_matrix_width, interpolatedCoordinates_width, annoNames_width, compress=FALSE);
	
	
	
	#####
	# Building average profile
	#####
	
	cat("Writing average profile\n");
	mean_values <- apply(result_matrix, MARGIN=2, mean);
	
	cairo_ps(filename=paste(output_folder, "profile-", expname,".ps",sep=""), width = 7, height = 7, bg = "transparent")
	plot(NULL, axes=F, main="G-quads", xlab="Position around center", ylab="G-quads occupancy", xlim=c(0,seq_length), ylim=c(0,max(mean_values)));
	lines(interpolatedCoordinates, mean_values, col="red");
	abline(v=seq_length/2,lty=2, col="green",lwd=2);
	axis(side=1,labels=c(-seq_length/2,-seq_length/4,"center", seq_length/4,seq_length/2), at=c(0,seq_length/4,seq_length/2,(seq_length/2)+(seq_length/4),seq_length), lwd=0, lwd.ticks=1, cex.axis=1.4);
	axis(side=2, lwd=0, lwd.ticks=1, cex.axis=1.4);
	box();
	dev.off();
}

################


##############
# MAIN
##############

getParams(paramsDefinition);

checkingOutputFolder(output_folder);

if(length(seq_length) != 1)
{
	stop("Sequence length should be a single numeric\n");
}

if(length(object_vec) != length(fasta_file_vec))
{
	stop("\n One object should be provided for each fasta file\n");
}

if(length(sorted_coord_gff_vec) != length(fasta_file_vec))
{
	stop("\n One gff file of coord should be given per fasta file\n");
}

for(i in 1:length(fasta_file_vec)) 
{
	#Reading fasta file
	fastaFile <- readDNAStringSet(fasta_file_vec[i], format="fasta", nrec=-1L, skip=0L, use.names=TRUE);
	
	#Reading object
	load(object_vec[i]);
	current_object <- genesSelected;
	rm(genesSelected);
	
	#Reading the sorted coordinates
	current_sorted_genes_gff <- read.table(sorted_coord_gff_vec[i], stringsAsFactors=FALSE);
	
	#Filtering the object regarding the ordered list of genes
	index_ordered <- genomicListgffFiltering(current_object, current_sorted_genes_gff, space_size_object, strandInformation = TRUE)
	current_object <- current_object[index_ordered];
	
	#Computing the mean value on a defined interval
	values_vec <- lapply(current_object, "[[", "probe.valueTSS");
	interval_filter <- c(intervalUpstream, intervalDownstream);
	values_vec <- filter_with_intervals(values_vec, interval_filter, binsize, profile_length, "TSS")
	mean_vec <- as.numeric(unlist(lapply(values_vec,mean)));
	
	cat("\nRetrieve the sequences\n");
	fastaFile <- as.list(as.character(fastaFile));
	
	cat("Looking at perfect G-quads in intrastrand fashion\n");
	
	coord_gquad_list_intrastrand_perfect <- mclapply(fastaFile, function(x){
				
				seq <- DNAString(x);
				pqs <- pqsfinder(seq,max_defects = 0)
				
				return(data.frame(start = start(pqs), end = start(pqs)+width(pqs), width=width(pqs), score=score(pqs)));
				
			}, mc.cores = nb_cpu);
	
	output_folder_intrastrand_perfect <- paste(output_folder, "intrastrand-perfect/", sep="");
	checkingOutputFolder(output_folder_intrastrand_perfect);
	exploring_gquads(coord_gquad_list_intrastrand_perfect, output_folder_intrastrand_perfect, expname, mean_vec, seq_length)
	
	rm(coord_gquad_list_intrastrand_perfect);
	gc();
	
	cat("Looking at G-quads in intrastrand fashion\n");
	
	coord_gquad_list_intrastrand <- mclapply(fastaFile, function(x){
				
				seq <- DNAString(x);
				pqs <- pqsfinder(seq)
				
				return(data.frame(start = start(pqs), end = start(pqs)+width(pqs), width=width(pqs), score=score(pqs)));
				
			}, mc.cores = nb_cpu);
	
	output_folder_intrastrand <- paste(output_folder, "intrastrand/", sep="");
	checkingOutputFolder(output_folder_intrastrand);
	exploring_gquads(coord_gquad_list_intrastrand, output_folder_intrastrand, expname, mean_vec, seq_length)
	
	rm(coord_gquad_list_intrastrand);
	gc();
	
	cat("Looking at perfect G-quads in interstrand fashion\n");
			
	coord_gquad_list_interstrand_perfect <- mclapply(fastaFile, function(x){
						
						seq <- DNAString(x);
						pqs <- pqsfinder(seq,max_defects = 0, strand = "+", use_default_scoring = FALSE, run_re = "G{3,6}|C{3,6}", custom_scoring_fn = isG4) 
						
						return(data.frame(start = start(pqs), end = start(pqs)+width(pqs), width=width(pqs), score=score(pqs)));
						
					}, mc.cores = nb_cpu);
	
	output_folder_interstrand_perfect <- paste(output_folder, "interstrand_perfect/", sep="");
	checkingOutputFolder(output_folder_interstrand_perfect);
	exploring_gquads(coord_gquad_list_interstrand_perfect, output_folder_interstrand_perfect, expname, mean_vec, seq_length);
	
	rm(coord_gquad_list_interstrand_perfect);
	gc();
	
	cat("Looking at G-quads in interstrand fashion\n");
	
	coord_gquad_list_interstrand <- mclapply(fastaFile, function(x){
				
				seq <- DNAString(x);
				pqs <- pqsfinder(seq, strand = "+", use_default_scoring = FALSE, run_re = "G{3,6}|C{3,6}", custom_scoring_fn = isG4) 
				
				return(data.frame(start = start(pqs), end = start(pqs)+width(pqs), width=width(pqs), score=score(pqs)));
				
			}, mc.cores = nb_cpu);
	
	output_folder_interstrand <- paste(output_folder, "interstrand/", sep="");
	checkingOutputFolder(output_folder_interstrand);
	exploring_gquads(coord_gquad_list_interstrand, output_folder_interstrand, expname, mean_vec, seq_length);
	
}

