##############################
# This script aims at developping a method for the ChIPSeqSpike package to 
# filter out reads that were aligned to both endogenous and exogenous genomes
# from BAM files.
# Descostes may 2018
##############################

source("class_and_methods.R")
source("common.R")

#########################
# PARAMETERS
#########################


infoFile_se <- "/home/descostes/Documents/test/test_chipseqspike/infopreprocess-single.csv"



#########################
# Methods and functions
#########################


.createBamIndex <- function(bam_file, verbose){
    if(!file.exists(paste0(bam_file, ".bai"))){
        if(verbose)
            message("\t\t Creating index for ", bam_file)
        indexBam(bam_file)
    }
}


.prepareOutFile <- function(bam_file){
    extension <- file_ext(bam_file)
    outFile <- paste0(file_path_sans_ext(bam_file), "-filtered.", extension)
    return(outFile)
}


.LoadFilterBamShortRead <- function(bam_file, out_file, rule, flg){
    
    filterBam(file = bam_file,
            destination = out_file,
            filter= rule,
            indexDestination=FALSE,
            ## Retrieve all BAM content (aligned)
            param=ScanBamParam(what=scanBamWhat(),
                    flag= flg))
}


.filterBAMAlignment <- function(object, endogenous_bam, exogenous_bam, 
        alignment_source, paired, verbose){
    
    .section(alignment_source, verbose, 1)
    
    ## Defining parameters to read bam files
    flagDef <- scanBamFlag(isUnmappedQuery = FALSE, isPaired= paired)
    pa <- ScanBamParam(flag = flagDef, what = "qname")
    
    if(!isTRUE(all.equal(length(endogenous_bam), 0))){
        if(!isTRUE(all.equal(length(exogenous_bam), 0))){
            
            bam_list <- list(endogenous_bam, exogenous_bam)
            
            ## Create index if does not exist
            invisible(lapply(bam_list, .createBamIndex, verbose))
            
            ## Preparing output files
            outfile_vec <- sapply(bam_list, .prepareOutFile)
            
            ## Reading BAM files
            if(verbose){
                message("\t Processing ", alignment_source, " reads")
                message("\t\t Reading endo and exo reads")
            }
            ID <- lapply(bam_list, function(x) scanBam(x, param=pa)[[1]][[1]])
            
            ## Retrieving common reads
            if(verbose)
                message("\t\t Retrieving common reads")
            commonID <- intersect(ID[[1]], ID[[2]])
            lcom <- length(commonID)
            if(verbose){
                lvec <- sapply(ID, length)
                perc <- sapply(lvec, function(lID){round((lcom*100)/lID)})
                message("\t\t\t", lcom, " reads in common (", perc[1], 
                        "% endo, ", perc[2], "% exo)")
                message("\t\t Filtering endo and exo BAM")
            }
            rm(ID)
            gc(verbose=FALSE)
            
            ## Filter BAM with filter rule
            filter_rule <- FilterRules(list(filterBAM = function(x)
                                !(x$qname %in% commonID)))
            ar <- list(filter_rule, flagDef)
            mapply(.LoadFilterBamShortRead, bam_list, outfile_vec, MoreArgs=ar)
            
            ## Update object fields
            if(isTRUE(all.equal(alignment_source, "Bowtie"))){
                endoBamBowtie(object) <- outfile_vec[1]
                exoBamBowtie(object) <- outfile_vec[2]
            }else{
                endoBamRsubread(object) <- outfile_vec[1]
                exoBamRsubread(object) <- outfile_vec[2]
            }
            return(object)
        }else{
            .endoNotExoError(object, alignment_source)
        }
    }
    
    return(object)
}


setMethod(
        
        f = "FilterCommonBam",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, paired = FALSE, verbose = TRUE){
            
            datasetList(theObject) <- lapply(getDatasetList(theObject), 
                    function(object){
                        
                        bamPrep <- .prepareBAM(object, verbose)
                        
                        ## Filtering reads
                        object <- .filterBAMAlignment(object, 
                                bamPrep[["enBamBow"]], bamPrep[["exBamBow"]],
                                "Bowtie", bamPrep[["paired"]], verbose)
                        object <- .filterBAMAlignment(object, 
                                bamPrep[["enBamRsub"]], bamPrep[["exBamRsub"]],
                                "Rsubread", bamPrep[["paired"]], verbose)
                        
                        return(object)
                    })
            
            return(theObject)
        })


setMethod(
        
        f = "AlignmentReport",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, paired = FALSE, verbose = TRUE){
            
            report <- lapply(getDatasetList(theObject), function(object){
                        
                        bam_vec <- c(getEndoBAMBowtie(object),
                                getExoBAMBowtie(object),
                                getEndoBAMRsubread(object),
                                getExoBAMRsubread(object))
                        
                        count_vec <- sapply(bam_vec, function(bam_file){
                                    if(verbose)
                                        message("Counting ", 
                                                basename(bam_file))
                                    pa <- ScanBamParam(
                                            scanBamFlag(isPaired= paired,
                                                    isUnmappedQuery = FALSE))
                                    result <- countBam(bam_file, 
                                            param=pa)$records
                                    return(result)
                                    
                                })
                        names(count_vec) <- basename(bam_vec)
                        return(count_vec)
                    })
            
            return(unlist(report))
        })


#########################
# MAIN
#########################


## Building the objects with complete alignments

csdsProcessSE <- spikePreprocessDataset(infoFile_se)

object1 <- getDatasetList(csdsProcessSE)[[1]]
endoBamBowtie(object1) <- ""
exoBamBowtie(object1) <- ""
object2 <- getDatasetList(csdsProcessSE)[[2]]
endoBamBowtie(object2) <- ""
exoBamBowtie(object2) <- ""
datasetList(csdsProcessSE) <- list(object1, object2)

# Filtering reads aligned to both genomes
FilterCommonBam(csdsProcessSE)
AlignmentReport(csdsProcessSE)
