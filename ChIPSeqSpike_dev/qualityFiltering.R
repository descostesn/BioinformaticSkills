##############################
# This script aims at developping a method for the ChIPSeqSpike package to 
# perform quality filtering on fastq files
# Descostes may 2018
##############################

source("class_and_methods.R")

#########################
# PARAMETERS
#########################

infoFile_pe <- ""
infoFile_se <- ""



#########################
# METHODS
#########################


setMethod(
        
        f = "fastqQualityReport",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, batch_size = 100000, kmer_length = 8, 
                outputFolder = NULL, width = 10, height = 10, verbose = FALSE){
            
            if(!(is.null(outputFolder) || file.exists(outputFolder)))
                dir.create(outputFolder, recursive = TRUE)
            
            fq_vec <- unlist(sapply(getDatasetList(theObject), 
                            function(object){
                                current_fastq1 <- getFASTQ1(object)
                                current_fastq2 <- getFASTQ2(object)
                                result <- current_fastq1
                                
                                if(!isTRUE(all.equal(length(current_fastq2), 0)))
                                    result <- c(result, current_fastq2)
                                
                                return(result)
                            }, simplify=FALSE))
            
            names(fq_vec) <- sapply(strsplit(basename(fq_vec), "\\."), "[",1)
            
            if(verbose)
                message("Computing stats on fastq files")
            fqlist <- seeFastq(fastq= fq_vec, batchsize= batch_size, 
                    klength= kmer_length)
            
            default_folder <- paste0(dirname(getFASTQ1(getDatasetList(
                                            theObject)[[1]])), "/")
            output_name <- paste0(
                    if(is.null(outputFolder)) default_folder else outputFolder,
                    "FastqQCreport.pdf")
            
            if(verbose)
                message("Writing ", output_name)
            pdf(file = output_name, width = width, height = height)
            seeFastqPlot(fqlist)
            dev.off()
        }
)


.filteringFastq <- function(fq, percentSeq, thresholdQual){
    
    qual_mat <- as(quality(fq), "matrix")
    nb_cycle <- (percentSeq * ncol(qual_mat))/100
    nb_cycle <- round(nb_cycle)
    validQual <- rowSums(qual_mat <= thresholdQual)
    fqFiltered <- fq[validQual <= nb_cycle]
    return(fqFiltered)
}


setMethod(
        
        f = "fastqQualityFiltering",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, chunk_size = 1000000, 
                percentSeq = 80, thresholdQual = 25, verbose = TRUE, ...){
            
            datasetList(theObject) <- lapply(getDatasetList(theObject), 
                    function(object){
                        fq1 <- getFASTQ1(object)
                        fq2 <- getFASTQ2(object)
                        paired <- !isTRUE(all.equal(length(fq2), 0))
                        extension <- paste(strsplit(basename(fq1), 
                                        "\\.")[[1]][-1], collapse=".")
                        output_file1 <- paste0(strsplit(fq1, "\\.")[[1]][1], 
                                "-filtered.", extension)
                        
                        if(paired){
                            b2 <- strsplit(fq2, "\\.")[[1]][1]
                            output_file2 <- paste0(b2, "-filtered.", extension)
                        }
                        
                        if(verbose){
                            message("Processing ", fq1, if(paired)
                                        paste0(" and ", fq2))
                        }
                        
                        totalNbReads <- qa(fq1)[["readCounts"]]$read
                        
                        if(totalNbReads < chunk_size)
                            chunk_size <- totalNbReads
                        
                        fq1stream <- FastqStreamer(fq1, chunk_size)
                        chunkCount <- 1
                        
                        if(paired){
                            totalNbReads2 <- qa(fq2)[["readCounts"]]$read
                            
                            if(!isTRUE(all.equal(totalNbReads, totalNbReads2)))
                                stop("The pair files do not have the same ",
                                        "number of reads")
                            fq2stream <- FastqStreamer(fq2, chunk_size)
                        }
                        
                        while(length(fq <- yield(fq1stream))) {
                            
                            fqFiltered <- .filteringFastq(fq, percentSeq, 
                                    thresholdQual)
                            idFiltered <- as.character(id(fqFiltered))
                            
                            if(paired){ 
                                fq2 <- yield(fq2stream)
                                fqFiltered2 <- .filteringFastq(fq2, percentSeq,
                                        thresholdQual)
                                idFiltered2 <- as.character(id(fqFiltered2))
                                
                                index1 <- idFiltered %in% idFiltered2
                                index2 <- idFiltered2 %in% idFiltered
                                
                                writeFastq(fqFiltered[index1], output_file1, 
                                        mode="a", ...)
                                writeFastq(fqFiltered[index2], output_file2, 
                                        mode="a", ...)
                            }else{
                                writeFastq(fqFiltered, output_file1, mode="a",
                                        ...)
                            }
                            
                            if(verbose){
                                amount <- (chunkCount*chunk_size)
                                percent <- round((amount*100)/totalNbReads)
                                cat("\r", percent, "%")
                                chunkCount <- chunkCount + 1
                            }
                        }
                        
                        close(fq1stream)
                        if(paired) close(fq2stream)
                        
                        totalNbReadsFilQ <- qa(output_file1)[["readCounts"]]
                        totalNbReadsFilQ <- totalNbReadsFilQ$read
                        
                        cat("\n\t", (totalNbReadsFilQ*100)/totalNbReads, 
                                "% of reads conserved after filtering\n")
                        
                        if(paired)
                            fastq2(object) <- output_file2
                        
                        fastq1(object) <- output_file1
                        
                        return(object)
                    })
            
            return(theObject)
        }
)



#########################
# MAIN
#########################


csdsProcessPE <- spikePreprocessDataset(infoFile_pe)
csdsProcessSE <- spikePreprocessDataset(infoFile_se)

fastqQualityReport(csdsProcessPE)
fastqQualityReport(csdsProcessSE, outputFolder = "./redirection/", 
        verbose = TRUE)

csdsProcessPE <- fastqQualityFiltering(csdsProcessPE, chunk_size = 1000000, 
        percentSeq = 80, thresholdQual = 25, verbose = TRUE)

csdsProcessSE <- fastqQualityFiltering(csdsProcessSE, chunk_size = 1000000, 
        percentSeq = 80, thresholdQual = 25, verbose = TRUE)

