##############################
# This script aims at developping a method for the ChIPSeqSpike package to 
# perform alignment with Rsubread and Rbowtie. The function should 
# enable to perform alignment with the 2 tools at the same time.
# Descostes ap 2018
##############################

source("class_and_methods.R")



#########################
# PARAMETERS
#########################

infoFile_pe <- ""
infoFile_se <- ""


#########################
# Methods
#########################


.buildIndexRsubread <- function(RsubreadIndex, genomeFasta, outputFolderIndex,
        endoExoString, verbose, ...){
    
    if(isTRUE(all.equal(length(RsubreadIndex), 0))){
        
        index_name <- file_path_sans_ext(basename(genomeFasta))
        if(!is.null(outputFolderIndex))
            index_name <- paste0(outputFolderIndex, index_name)
        
        if(verbose)
            message("\t Building ", endoExoString, " genome index with ",
                    "Rsubread")
        buildindex(basename = index_name, reference = genomeFasta, ...)
        RsubreadIndex <- if(!is.null(outputFolderIndex)) index_name
                else paste0(dirname(genomeFasta), "/", index_name)
        
        return(RsubreadIndex)
    }else{
        return(character())
    }
}


.RsubreadAlignment <- function(readFile1, readFile2, endoExoString, 
        RsubreadIndex, input_format, outputFolder, phred64, unique_mapping, 
        nbCpu, genomeName, verbose, ...){
    
    if(verbose)
        if(!is.null(readFile2))
            message("\t Performing paired-end ", endoExoString, 
                    " Rsubread alignment of ", readFile1, " on ", genomeName)
        else
            message("\t Performing ", endoExoString, " alignment of ", 
                    readFile1)
    
    result_file <- paste0(if(!is.null(outputFolder)) outputFolder else
                        paste0(dirname(readFile1), "/"),
            strsplit(basename(readFile1),"\\.")[[1]][1], "_",
            genomeName, if(unique_mapping) "_unique" else 
                        "_multiple", 
            if(isTRUE(all.equal(length(readFile2), 0))) "_SE" else
                        "_PE", "_subread.BAM")
    
    align(index = RsubreadIndex,
            readfile1 = readFile1,
            readfile2 = if(!isTRUE(all.equal(length(readFile2), 0))) readFile2 
                    else NULL,
            type="dna",
            input_format = input_format,
            output_format = "BAM",
            output_file = result_file,
            # offset value added to Phred quality scores of read bases
            phredOffset = if(phred64) 64 else 33,
            unique = unique_mapping,
            nthreads = nbCpu, ...)
    
    return(result_file)
}


.BowtieAlignment <- function(sampleFile_path, genomeFasta, readFile1, 
        readFile2, outputFolder, nbCpu, genomeName, endoExoString, 
        unique_mapping, multireadThreshold, verbose, ...){
    
    if(isTRUE(all.equal(length(readFile2), 0))){ ## Single-end
        
        sampleFile <- paste("FileName", "SampleName", sep="\t")
        sampleFile <- c(sampleFile, paste(readFile1, 
                        strsplit(basename(readFile1),
                                "\\.")[[1]][1], sep="\t"))
    }else{ ## Paired-end
        paired_name <- paste0(strsplit(basename(readFile1),
                        "\\.")[[1]][1], 
                strsplit(basename(readFile2),"\\.")[[1]][1])
        sampleFile <- paste("FileName1", "FileName2", "SampleName",
                sep="\t")
        sampleFile <- c(sampleFile, paste(readFile1,
                        readFile2,
                        paired_name, sep="\t"))
    }
    
    sampleFile_path <- paste0(dirname(readFile1), 
            "/tmp-sample.txt")
    write(sampleFile, file= sampleFile_path)
    
    if(verbose)
        if(!is.null(readFile2))
            message("\t Performing paired-end ", endoExoString, 
                    " Bowtie alignment of ", readFile1, " on ", genomeName)
        else
            message("\t Performing ", endoExoString, " alignment of ", 
                    readFile1)
    
    
    result_bowtie1 <- qAlign(sampleFile = sampleFile_path,
            genome = genomeFasta,
            maxHits= if(unique_mapping) 1 else multireadThreshold, 
            paired= if(isTRUE(all.equal(length(readFile2), 0))) "no" else "fr",
            alignmentsDir= if(!is.null(outputFolder)) outputFolder else NULL, 
            clObj= makeCluster(nbCpu), ...)
    
    # Organizing and renaming files
    unlink(sampleFile_path)
    unlink("QuasR*")
    current_dir <- getwd()
    on.exit(setwd(current_dir))
    if(!is.null(outputFolder)) 
        setwd(outputFolder) 
    else 
        setwd(dirname(readFile1))
    tmp_txt <- list.files("./", ".bam.txt")
    tmp_bam <- paste(strsplit(tmp_txt, "\\.")[[1]][1:2], collapse=".")
    unlink(tmp_txt)
    tmp_bai <- paste0(tmp_bam, ".bai")
    md5nb <- strsplit(strsplit(tmp_bam, "_")[[1]][2], "\\.")[[1]][1]
    
    new_name <- paste0(genomeName, if(unique_mapping) "_unique" else 
                        "_multiple", 
            if(isTRUE(all.equal(length(readFile2), 0))) "_SE" else
                        "_PE", "_bowtie")
    
    file.rename(tmp_bam, sub(md5nb, new_name, tmp_bam))
    file.rename(tmp_bai, sub(md5nb, new_name, tmp_bai))
    result_file <- paste0(if(!is.null(outputFolder)) outputFolder else 
                        paste0(dirname(readFile1), "/"),
            sub(md5nb, new_name, tmp_bam))
    rm(list=c("tmp_txt", "tmp_bam", "tmp_bai", "md5nb"))
    
    return(result_file)
}


setMethod(
        
        f = "ChIPReadsAlign",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, unique_mapping = TRUE, 
                multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
                outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
                Bowtie = FALSE, verbose = FALSE, ...){
            
            datasetList(theObject) <- lapply(getDatasetList(theObject), 
                    function(object){
                        ChIPReadsAlign(object, unique_mapping, 
                                multireadThreshold, nbCpu, phred64, 
                                outputFolderIndex, outputFolder, Rsubread,
                                Bowtie, verbose, ...)
                    })
            
            return(theObject)
        }
)


setMethod(
        
        f = "ChIPReadsAlign",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, unique_mapping = TRUE, 
                multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
                outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
                Bowtie = FALSE, verbose = FALSE, ...){
            
            if(!(is.null(outputFolderIndex) || file.exists(outputFolderIndex)))
                dir.create(outputFolderIndex, recursive = TRUE)
            
            if(!(is.null(outputFolder) || file.exists(outputFolder)))
                dir.create(outputFolder, recursive = TRUE)
            
            
            endoGenomeFasta <- getFasta(theObject)
            exoGenomeFasta <- getExoFasta(theObject)
            current_fastq1 <- getFASTQ1(theObject)
            current_fastq2 <- getFASTQ2(theObject)
            endoGenome <- getGenome(theObject)
            exoGenome <- getExoGenome(theObject)
            
            
            if(Rsubread){
                if(verbose)
                    message("Performing alignment with Rsubread")
                
                endoRsubreadIndex <- getRsubreadIndex(theObject)
                exoRsubreadIndex <- getExoRsubreadIndex(theObject)
                
                ## Building index if necessary
                if(isTRUE(all.equal(length(endoRsubreadIndex), 0)) ||
                        isTRUE(all.equal(length(exoRsubreadIndex), 0))){
                    
                    if(is.null(outputFolderIndex)){
                        currentwd <- getwd()
                        setwd(dirname(endoGenomeFasta))
                    }
                    indexRsubread(theObject) <- .buildIndexRsubread(
                            endoRsubreadIndex, endoGenomeFasta,
                            outputFolderIndex, "endogenous", verbose, ...)
                    if(is.null(outputFolderIndex))
                        setwd(dirname(exoGenomeFasta))
                    exoIndexRsubread(theObject) <- .buildIndexRsubread(
                            exoRsubreadIndex, exoGenomeFasta,
                            outputFolderIndex, "exogenous", verbose, ...)
                    
                    if(is.null(outputFolderIndex))
                        setwd(currentwd)
                }
                
                endoRsubreadIndex <- getRsubreadIndex(theObject)
                exoRsubreadIndex <- getExoRsubreadIndex(theObject)
                
                ## Determining input format
                extension <- paste(strsplit(basename(current_fastq1),
                                "\\.")[[1]][-1], collapse=".")
                if(identical(extension, "fq") || identical(extension, "fastq"))
                    input_format <- "FASTQ"
                else
                    input_format <- "gzFASTQ"
                
                ## Performing alignment with Rsubread on both endogenous and 
                ## exogenous genomes
                endoBamRsubread(theObject) <- character()
                endoBamRsubread(theObject) <- .RsubreadAlignment(
                        current_fastq1, current_fastq2, "endogenous", 
                        endoRsubreadIndex, input_format, outputFolder, phred64,
                        unique_mapping, nbCpu, endoGenome, verbose, ...)
                
                exoBamRsubread(theObject) <- character()
                exoBamRsubread(theObject) <- .RsubreadAlignment(current_fastq1,
                        current_fastq2, "exogenous", exoRsubreadIndex, 
                        input_format, outputFolder, phred64, unique_mapping, 
                        nbCpu, exoGenome, verbose, ...)
            }
            
            if(Bowtie){
                
                ## Performing alignment with Bowtie on both endogenous and 
                ## exogenous genomes
                endoBamBowtie(theObject) <- .BowtieAlignment(sampleFile_path, 
                        endoGenomeFasta, current_fastq1, current_fastq2, 
                        outputFolder, nbCpu, endoGenome, "endogenous", 
                        unique_mapping, multireadThreshold, verbose, ...)
                
                exoBamBowtie(theObject) <- .BowtieAlignment(sampleFile_path, 
                        exoGenomeFasta, current_fastq1, current_fastq2, 
                        outputFolder, nbCpu, exoGenome, "exogenous", 
                        unique_mapping, multireadThreshold, verbose, ...)
            }
            
            return(theObject)
        }
)




#########################
# MAIN
#########################


csdsProcessPE <- spikePreprocessDataset(infoFile_pe)
csdsProcessSE <- spikePreprocessDataset(infoFile_se)


## testing single-end - current folder - unireads
csdsProcessSECFUR <- ChIPReadsAlign(csdsProcessSE, unique_mapping = TRUE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

## testing single-end - current folder - multireads
csdsProcessSECFMR <- ChIPReadsAlign(csdsProcessSE, unique_mapping = FALSE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)


## testing paired-end - current folder - unireads
csdsProcessPECFUR <- ChIPReadsAlign(csdsProcessPE, unique_mapping = TRUE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

## testing paired-end - current folder - multireads
csdsProcessPECFMR <- ChIPReadsAlign(csdsProcessPE, unique_mapping = FALSE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)



## testing single-end - other output folder - unireads
csdsProcessSEOFUR <- ChIPReadsAlign(csdsProcessSE, unique_mapping = TRUE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = "/home/descostes/Documents/test/test_chipseqspike/redirection/SEUN/", Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

## testing single-end - other folder - multireads
csdsProcessSEOFMR <- ChIPReadsAlign(csdsProcessSE, unique_mapping = FALSE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = "/home/descostes/Documents/test/test_chipseqspike/redirection/SEMR/", Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

## testing paired-end - other output folder - unireads
csdsProcessPEOFUR <- ChIPReadsAlign(csdsProcessPE, unique_mapping = TRUE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = "/home/descostes/Documents/test/test_chipseqspike/redirection/PEUN/", Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

## testing paired-end - other folder - multireads
csdsProcessPEOFMR <- ChIPReadsAlign(csdsProcessPE, unique_mapping = FALSE, 
        multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
        outputFolderIndex = NULL, outputFolder = "/home/descostes/Documents/test/test_chipseqspike/redirection/PEMR/", Rsubread = TRUE,
        Bowtie = TRUE, verbose = TRUE)

