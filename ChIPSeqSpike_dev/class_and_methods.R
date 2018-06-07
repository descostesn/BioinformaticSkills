
#library("Rsubread")
library("QuasR")
library("parallel")
library("tools")
library("systemPipeR")
library("ShortRead")
library("Rsamtools")
library("GenomicAlignments")


######################################################################
## !!!!!!!!!!!!! ALREADY IN THE PACKAGE, DO NOT COPY BACK
######################################################################


.validatePath <- function(value){
    
    if(!file.exists(value))
        stop(value, " is not a valid path")
    
    if(length(gregexpr(pattern="\\.", basename(value))[[1]]) > 1)
        stop(value, " should contain only one point for the extension.")
    
}

.readInfoFile <- function(infoFile){
    
    if(isTRUE(all.equal(file_ext(infoFile), "csv"))){
        info_table <- read.csv(infoFile, stringsAsFactors = FALSE)
    }else if(isTRUE(all.equal(file_ext(infoFile), "txt"))){
        info_table <- read.table(infoFile, header = TRUE, sep = "\t", 
                stringsAsFactors = FALSE)
    }else
        stop("The info file should be in csv or txt format.")
    
    return(info_table)
}

.checkInfoFile <- function(info_table, expected_nb_col, colnames_vec){
    
    if(!(isTRUE(all.equal(ncol(info_table), expected_nb_col)) && 
                isTRUE(all.equal(colnames(info_table), 
                                colnames_vec))))
        stop("Columns are missing in your info file or one column has the ",
                "wrong name. Your file should contain the following column ",
                "names: ", paste(colnames_vec, collapse=","))
    
    if(isTRUE(all.equal(nrow(info_table), 0)))
        stop("File is empty.")
}

setGeneric(
        
        name = "datasetList<-",
        
        def = function(theObject, value){
            
            standardGeneric("datasetList<-")
        },
        signature = "theObject")

##################################################################

#########################
# CLASSES
#########################


ChIPSeqSpikeProcessList <- setClass(
        
        Class = "ChIPSeqSpikeProcessList",
        
        slots = c(alignList = "list"),
        
        validity = function(object){
            
            if(length(object@alignList) < 1)
                stop("must contain at least one ChIPSeqSpikeProcess object.")
            else if(!all(vapply(object@alignList, function(x) 
                                inherits(x, "ChIPSeqSpikeProcess"), logical(1))
            ))
                stop("all objects in list must be of ChIPSeqSpikeProcess", 
                        " type.")
        })


ChIPSeqSpikeProcess <- setClass(
        
        Class = "ChIPSeqSpikeProcess",
        
        slots = c(experimentFile_reads1 = "character",
                experimentFile_reads2 = "character",
                endogenousGenomeFile = "character",
                endogenousGenomeName = "character",
                endogenousRsubreadIndex = "character",
                exogenousGenomeFile = "character",
                exogenousGenomeName = "character",
                exogenousRsubreadIndex = "character",
                endoBamFileSubread = "character",
                endoBamFileBowtie = "character",
                exoBamFileSubread = "character",
                exoBamFileBowtie = "character"),
        
        validity = function(object){
            
            .validateGenome_or_FastqFile(object@experimentFile_reads1, "fastq")
            .validateGenome_or_FastqFile(object@endogenousGenomeFile, "fasta")
            .validateGenome_or_FastqFile(object@exogenousGenomeFile, "fasta")
            
            if(!isTRUE(all.equal(length(object@endoBamFileSubread),0)))
                .validatePath(object@endoBamFileSubread)
            
            if(!isTRUE(all.equal(length(object@exoBamFileSubread),0)))
                .validatePath(object@exoBamFileSubread)
            
            if(!isTRUE(all.equal(length(object@endoBamFileBowtie),0)))
                .validatePath(object@endoBamFileBowtie)
            
            if(!isTRUE(all.equal(length(object@exoBamFileBowtie),0)))
                .validatePath(object@exoBamFileBowtie)
            
            if(!isTRUE(all.equal(length(object@endogenousRsubreadIndex),0)))
                .validateRSubReadIndex(object@endogenousRsubreadIndex)
            
            if(!isTRUE(all.equal(length(object@exogenousRsubreadIndex), 0)))
                .validateRSubReadIndex(object@exogenousRsubreadIndex)
        })


.validateGenome_or_FastqFile <- function(value, format){
    
    if(!file.exists(value))
        stop(value, " is not a valid path. fasta and fastq are required.")
    
    nb_point <- length(gregexpr(pattern="\\.", value)[[1]])
    extension <- strsplit(basename(value), "\\.")[[1]][2]
    no_underscore_path <- strsplit(dirname(value), "_")
    
    if(nb_point > 2 || !(identical(extension, if(identical(format, "fasta")) 
                                    "fa" else "fq") || 
                identical(extension, if(identical(format, "fasta")) "fasta"
                                else "fastq"))){
        if(identical(format, "fasta"))
            stop("Fasta files accepted extension are .fa or ", 
                    ".fasta. File path cannot contain points.")
        else
            stop("Fastq files accepted extension are .fq or ", 
                    ".fastq. File path cannot contain points.")
    }
}

.validateRSubReadIndex <- function(value){
    
    if(!(file.exists(paste0(value, ".00.b.array")) &&
                file.exists(paste0(value, ".00.b.tab")) &&
                file.exists(paste0(value, ".files")) &&
                file.exists(paste0(value, ".log")) &&
                file.exists(paste0(value, ".reads"))
                ))
        stop(value, " should contain the following files: ",
                paste0(value, ".00.b.array, "),
                paste0(value, ".00.b.tab, "),
                paste0(value, ".files, "),
                paste0(value, ".log, and "),
                paste0(value, ".reads"))
}


#########################
# GENERICS DECLARATION
#########################


setGeneric(
        
        name = "getDatasetList",
        
        def = function(theObject){
            
            standardGeneric("getDatasetList")
        },
        signature = "theObject")

setGeneric(
        
        name = "getFasta",
        
        def = function(theObject){
            
            standardGeneric("getFasta")
        },
        signature = "theObject")

setGeneric(
        
        name = "getExoFasta",
        
        def = function(theObject){
            
            standardGeneric("getExoFasta")
        },
        signature = "theObject")

setGeneric(
        
        name = "getRsubreadIndex",
        
        def = function(theObject){
            
            standardGeneric("getRsubreadIndex")
        },
        signature = "theObject")

setGeneric(
        
        name = "getExoRsubreadIndex",
        
        def = function(theObject){
            
            standardGeneric("getExoRsubreadIndex")
        },
        signature = "theObject")

setGeneric(
        
        name = "getFASTQ1",
        
        def = function(theObject){
            
            standardGeneric("getFASTQ1")
        },
        signature = "theObject")

setGeneric(
        
        name = "getFASTQ2",
        
        def = function(theObject){
            
            standardGeneric("getFASTQ2")
        },
        signature = "theObject")

setGeneric(
        
        name = "getGenome",
        
        def = function(theObject){
            
            standardGeneric("getGenome")
        },
        signature = "theObject")

setGeneric(
        
        name = "getExoGenome",
        
        def = function(theObject){
            
            standardGeneric("getExoGenome")
        },
        signature = "theObject")

setGeneric(
        
        name = "getExoBAMRsubread",
        
        def = function(theObject){
            
            standardGeneric("getExoBAMRsubread")
        },
        signature = "theObject")

setGeneric(
        
        name = "getEndoBAMRsubread",
        
        def = function(theObject){
            
            standardGeneric("getEndoBAMRsubread")
        },
        signature = "theObject")

setGeneric(
        
        name = "getExoBAMBowtie",
        
        def = function(theObject){
            
            standardGeneric("getExoBAMBowtie")
        },
        signature = "theObject")

setGeneric(
        
        name = "getEndoBAMBowtie",
        
        def = function(theObject){
            
            standardGeneric("getEndoBAMBowtie")
        },
        signature = "theObject")

#setReplaceMethod


setGeneric(
        
        name = "indexRsubread<-",
        
        def = function(theObject, value){
            
            standardGeneric("indexRsubread<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "exoIndexRsubread<-",
        
        def = function(theObject, value){
            
            standardGeneric("exoIndexRsubread<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "endoBamRsubread<-",
        
        def = function(theObject, value){
            
            standardGeneric("endoBamRsubread<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "exoBamRsubread<-",
        
        def = function(theObject, value){
            
            standardGeneric("exoBamRsubread<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "endoBamBowtie<-",
        
        def = function(theObject, value){
            
            standardGeneric("endoBamBowtie<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "exoBamBowtie<-",
        
        def = function(theObject, value){
            
            standardGeneric("exoBamBowtie<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "fastq1<-",
        
        def = function(theObject, value){
            
            standardGeneric("fastq1<-")
        },
        signature = "theObject")

setGeneric(
        
        name = "fastq2<-",
        
        def = function(theObject, value){
            
            standardGeneric("fastq2<-")
        },
        signature = "theObject")

# Alignment method

setGeneric(
        
        name = "ChIPReadsAlign",
        
        def = function(theObject, unique_mapping = TRUE, 
                multireadThreshold = 200, nbCpu = 1, phred64 = FALSE, 
                outputFolderIndex = NULL, outputFolder = NULL, Rsubread = TRUE,
                Bowtie = FALSE, verbose = FALSE, ...){
            
            standardGeneric("ChIPReadsAlign")
        },
        signature = "theObject")

# Filtering methods


setGeneric(
        
        name = "fastqQualityReport",
        
        def = function(theObject, batch_size = 100000, kmer_length = 8, 
                outputFolder = NULL, width = 10, height = 10, verbose = FALSE){
            
            standardGeneric("fastqQualityReport")
        },
        signature = "theObject")


setGeneric(
        
        name = "fastqQualityFiltering",
        
        def = function(theObject, chunk_size = 1000000, percentSeq = 80, 
                thresholdQual = 25, verbose = TRUE, ...){
            
            standardGeneric("fastqQualityFiltering")
        },
        signature = "theObject")

# Common BAM and alignment report


setGeneric(
        
        name = "FilterCommonBam",
        
        def = function(theObject, paired = FALSE, verbose = TRUE){
            
            standardGeneric("FilterCommonBam")
        },
        signature = "theObject")


setGeneric(
        
        name = "AlignmentReport",
        
        def = function(theObject, paired = FALSE, verbose = TRUE){
            
            standardGeneric("AlignmentReport")
        },
        signature = "theObject")


setGeneric(
        
        name = "filterDuplicates",
        
        def = function(theObject, aligners = c("Bowtie", "Rsubread"), 
                nbCPU = 1, reportPDF = TRUE, outputFolder = NULL, num_sim = 1,
                dupThreshold = NULL, blackListEndo = NULL, blackListExo = NULL,
                verbose = TRUE, ...){
            
            standardGeneric("filterDuplicates")
        },
        signature = "theObject")


#########################
# ACCESSORS
#########################

setMethod(
        
        f = "getDatasetList",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject){
            
            return(theObject@alignList)
        })

setMethod(
        
        f = "getFasta",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@endogenousGenomeFile)
        })

setMethod(
        
        f = "getExoFasta",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@exogenousGenomeFile)
        })

setMethod(
        
        f = "getRsubreadIndex",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@endogenousRsubreadIndex)
        })

setMethod(
        
        f = "getExoRsubreadIndex",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@exogenousRsubreadIndex)
        })

setMethod(
        
        f = "getFASTQ1",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@experimentFile_reads1)
        })

setMethod(
        
        f = "getFASTQ2",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@experimentFile_reads2)
        })

setMethod(
        
        f = "getGenome",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@endogenousGenomeName)
        })

setMethod(
        
        f = "getExoGenome",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@exogenousGenomeName)
        })

setMethod(
        
        f = "getExoBAMRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@exoBamFileSubread)
        })

setMethod(
        
        f = "getEndoBAMRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@endoBamFileSubread)
        })

setMethod(
        
        f = "getExoBAMBowtie",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@exoBamFileBowtie)
        })

setMethod(
        
        f = "getEndoBAMBowtie",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject){
            
            return(theObject@endoBamFileBowtie)
        })


#########################
# REPLACE METHODS
#########################

setReplaceMethod(
        
        f = "datasetList",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, value){
            
            theObject@alignList <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "indexRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@endogenousRsubreadIndex <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "exoIndexRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@exogenousRsubreadIndex <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "endoBamRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@endoBamFileSubread <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "exoBamRsubread",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@exoBamFileSubread <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "endoBamBowtie",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@endoBamFileBowtie <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "exoBamBowtie",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@exoBamFileBowtie <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "fastq1",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@experimentFile_reads1 <- value
            validObject(theObject)
            return(theObject)
        })

setReplaceMethod(
        
        f = "fastq2",
        
        signature = "ChIPSeqSpikeProcess",
        
        definition = function(theObject, value){
            
            theObject@experimentFile_reads2 <- value
            validObject(theObject)
            return(theObject)
        })


#########################
# CONSTRUCTORS
#########################



spikePreprocessDataset <- function(infoFile){
    
    info_table <- .readInfoFile(infoFile)
    colnames_vec <- c("file1", "file2", "endoFasta", "endoName", 
            "endoRsubreadIndex", "exoFasta", "exoName", "exoRsubreadIndex")
    .checkInfoFile(info_table, 8, colnames_vec)
    
    list_dataset <- split(info_table, factor(info_table$file1, 
                    levels = info_table$file1))
    names(list_dataset) <- paste0("dataset", 
            seq_len(length(list_dataset)))
    
    ChIPSeqSpikeProcessList(list_dataset)
}


ChIPSeqSpikeProcessList <- function(dataset_list){
    
    ChIPSeqSpikeProcessCollection <- lapply(dataset_list, 
            function(dataset_table){
                
                endoRsubreadIndex <- dataset_table[, "endoRsubreadIndex"]
                exoRsubreadIndex <- dataset_table[, "exoRsubreadIndex"]
                
                return(
                        new("ChIPSeqSpikeProcess",
                                experimentFile_reads1 = 
                                        dataset_table[, "file1"],
                                experimentFile_reads2 = 
                                        if(!is.na(dataset_table[, "file2"]))
                                            dataset_table[, "file2"] else
                                            character(),
                                endogenousGenomeFile = 
                                        dataset_table[, "endoFasta"],
                                endogenousGenomeName = 
                                        dataset_table[, "endoName"],
                                endogenousRsubreadIndex = 
                                        if(is.na(endoRsubreadIndex)) 
                                            character() else endoRsubreadIndex,
                                exogenousGenomeFile = 
                                        dataset_table[, "exoFasta"],
                                exogenousGenomeName = 
                                        dataset_table[, "exoName"],
                                exogenousRsubreadIndex = 
                                        if(is.na(exoRsubreadIndex)) 
                                            character() else exoRsubreadIndex,
                                endoBamFileSubread = character(),
                                endoBamFileBowtie = character(),
                                exoBamFileSubread = character(),
                                exoBamFileBowtie = character())
                )
            })
    
    new("ChIPSeqSpikeProcessList", 
            alignList = ChIPSeqSpikeProcessCollection)
}


