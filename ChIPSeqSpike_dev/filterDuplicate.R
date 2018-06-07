##############################
# This script aims at developping a method for the ChIPSeqSpike package to 
# filter duplicates of PCR (with htseqtools) and removing reads found into 
# blacklists (with ??).
# Descostes may 2018
##############################

#library("htSeqTools"), is going to be removed in the next release of bioconductor

library("gridExtra")

source("class_and_methods.R")
source("common.R")
source("htseqtoolsmethods.R")
source("fdrEnrichedCounts_htSeqTools.R")

#########################
# PARAMETERS
#########################

infoFile_pe <- ""
infoFile_se <- ""


nbCPU <- 1

#mock

#########################
# FUNCTION
#########################


## ---------- Subfunctions associated with .filterPCRduplicates

.computeGenomicAlignments <- function(endoExolist, pairedEnd, verb){
    
    ga <- lapply(endoExolist, function(bam_files){
                return(lapply(bam_files, function(bamfile){
                                    if(verb)
                                        message("Reading and converting ", 
                                                bamfile)
                                    if(pairedEnd) 
                                        granges(readGAligmentPairs(bamfile))
                                    else granges(readGAlignments(bamfile))
                                }))
            })
    ga <- lapply(ga, GRangesList)
    return(ga)
}

.SsdAndGini <- function(gaObject, msg, verbose, f, ...){
    
    if(verbose)
        message("\t", msg)
    
    result <- lapply(gaObject, function(ga_element){
                lapply(ga_element, function(x) f(x,...))
            })
    return(result)
}

.computeExpNames <- function(endoExolist){
    
    names_vec <- vapply(endoExolist, function(bams){names(bams)}, 
            character(2), USE.NAMES = FALSE)
    return(as.vector(names_vec))
}

.computeSSD <- function(gaObj, namesv, verb){
    
    ssd_vec <- .SsdAndGini(gaObj, "Computing ssd", verb, 
            ssdCoverage)
    ssd_vec <- round(unlist(ssd_vec), digits = 2)
    names(ssd_vec) <- namesv
    return(ssd_vec)
}

.computeNRepeats <- function(gaObj, verb){
    
    if(verb)
        message("\t Determining numbers of duplicates")
    
    nrepeats_list <- lapply(gaObj, function(bamfiles){
                duplist <- lapply(bamfiles, function(bam){
                            return(tabDuplReads(bam, nbCPU))})
                return(duplist)
            })
    return(nrepeats_list)
}

.computeGini <- function(gaObj, ncpu, namesv, verb, sim){
    
    ## The returned list contains a list for each bam file:
    ## Each element of the list has a list with gini coeff as first 
    ## element and a list of two functions as second element.
    ## The two functions are without arguments and enable to plot: 
    ## 1) proportions  and 2) the lorenz curve
    gini_list <- .SsdAndGini(gaObj, "Computing Gini Coverage", 
            verb, giniCoverage, mc.cores = ncpu, numSim = sim)
    
    # Building gini index table
    # First vapply separate endogenous and exogenous
    # Second vapply separates experiments
    gini_table <- lapply(gini_list, function(genomeType){
                
                typeTable <- vapply(genomeType, function(experiment){
                            return(experiment[["ginicoeff"]])
                        }, numeric(2))
                return(t(typeTable))
            })
    gini_table <- do.call(rbind, gini_table)
    
    # Building proportion function vector
    gini_propFUN <- lapply(gini_list, function(genomeType){
                lapply(genomeType, function(experiment){
                            return(experiment[["giniplots"]][[1]])
                        })
            })
    
    # Building lorenz curve function vector
    gini_lorenzFUN <- lapply(gini_list, function(genomeType){
                lapply(genomeType, function(experiment){
                            return(experiment[["giniplots"]][[1]])
                        })
            })
    
    return(list(table = gini_table, proportion = gini_propFUN, 
                    lorenzCurve = gini_lorenzFUN))
}

.convertToIRanges <- function(experiment){
    return(split(ranges(experiment), 
                    factor(as.vector(seqnames(experiment)), 
                            levels = seqlevels(experiment))))
}

.convertToIRangesList <- function(gaObject){
    
    result <- lapply(gaObject, function(genomeType){
                
                ## IRL = IRangesList
                genomeTypeIRL <- lapply(genomeType, function(experiment){
                            IRanges_list <- .convertToIRanges(experiment)
                            return(IRanges_list)
                        })
                return(genomeTypeIRL)
            })
    return(result)
}

.filterDuplReads <- function(gaIRanges, countrepeatslist, maxRepeats, nbCPU, 
        verbose, ...){
    
    filteredgaIRanges <- mapply(function(genomeType, countrepeatsElement){
                
                filteredType <- mapply(function(experiment, countrepeats){
                            filterDuplReads(experiment, countrepeats, 
                                    maxRepeats, nbCPU, verbose, ...)
                        }, genomeType, countrepeatsElement, SIMPLIFY = FALSE)
                return(filteredType)
            }, gaIRanges, countrepeatslist, SIMPLIFY = FALSE)
    
    return(filteredgaIRanges)
}

.retrieveFilteredBam <- function(object){
    result <- lapply(object, function(genomeType){
                result2 <- lapply(genomeType, function(experiment){
                            return(experiment[["filteredobject"]])
                        })
                return(result2)
            })
    return(result)
}

## .thresholdFUN and .totalReadsFUN are used as arguments of 
## .retrieveFilteringInfo 
.thresholdFUN <- function(x) return(x[["threshold"]])
.totalReadsFUN <- function(x) return(sum(sapply(x, length)))
.retrieveFilteringInfo <- function(object, f){
    
    result <- unlist(lapply(object, function(genomeType){
                        result2 <- sapply(genomeType, 
                                function(experiment){
                                    return(f(experiment))
                                }
                        )
                        return(result2)
                    }
            ))
    return(result)
}

.testTresholdEffect <- function(gaIRanges, estimated_thresholds, 
        verbose){
    
    thresholdTestVec <- seq_len(max(estimated_thresholds))
    
    resultNbReadsList <- lapply(gaIRanges, function(genomeType){
                
                nameExps <- names(genomeType)
                
                expsNbReads <- mapply(function(experiment, 
                                nameExp){
                            
                            if(verbose)
                                message("\t", nameExp)
                            
                            result <- filterDuplReads(experiment, 
                                    nrepeats_list, dupThreshold, 
                                    nbCPU, 
                                    testThresholdVec = 
                                            thresholdTestVec, 
                                    verbose)
                            readsKeptVec <- sapply(result, 
                                    .totalReadsFUN)
                            names(readsKeptVec) <- thresholdTestVec
                            return(readsKeptVec)
                        }, genomeType, nameExps, SIMPLIFY = FALSE)
                return(expsNbReads)
            })
    
    return(resultNbReadsList)
}

.filterBlackList <- function(gaIRanges_filtered, blackLists_list, 
        verbose){
    
    result <- mapply(function(genomeType, blackList, name, pairedEnd){
                
                if(!is.null(blackList)){
                    if(verbose)
                        message("Reading ", name, " black list")
                    formatAnno <- strsplit(basename(blackList), "\\.")[[1]][2]
                    blacklistGRanges <- import(blackList, format = formatAnno)
                    blackListIRanges <- .convertToIRanges(blacklistGRanges)
                    
                    filteredList <- lapply(genomeType, function(experiment){
                                
                                experiment <- IRangesList(experiment)
                                resOv <- findOverlaps(experiment, 
                                        blackListIRanges)
                                filteredExp <- mapply(function(expchr, ovchr){
                                            indRemov <- queryHits(ovchr)
                                            if(!isTRUE(
                                                    all.equal(length(indRemov),
                                                            0)
                                            ))
                                                return(expchr[-indRemov, ])
                                            else return(expchr)
                                        }, experiment, resOv)
                            })
                    return(filteredList)
                }else{
                    return(genomeType)
                }
                
            }, gaIRanges_filtered, blackLists_list, names(blackLists_list),
            SIMPLIFY = FALSE)
    return(result)
}

## ----------


.filterPCRduplicates <- function(object, endoExoBamlist, paired, nbCPU, 
        reportPDF, outputFolder, num_sim, dupThreshold, blackListEndo, 
        blackListExo, verbose, ...){
    
    .section("Filtering PCR duplicates", verbose, 2)
    
    if(is.null(outputFolder)) outputFolder <- dirname(getFASTQ1(object))
    
    ga <- .computeGenomicAlignments(endoExoBamlist, paired, verbose)
    gaIRanges <- .convertToIRangesList(ga)
    
    nrepeats_list <- .computeNRepeats(gaIRanges, verbose)
    
    if(verbose)
        message("Filtering duplicates")
    result_filtering <- .filterDuplReads(gaIRanges, nrepeats_list, 
            dupThreshold, nbCPU, verbose, ...)
    
    if(verbose)
        message("Retrieving filtered objects")
    gaIRanges_filtered <- .retrieveFilteredBam(result_filtering)
    
    browser()
    do_blacklist <- !is.null(blackListEndo) || !is.null(blackListExo)
    
    if(do_blacklist){
        
        if(verbose)
            .section("Filtering out reads found in ENCODE blacklists", 
                    verbose, 2)
        
        blackList_list <- list(endo = blackListEndo, exo = blackListExo)
        gaIRanges_filtered_blackList <- .filterBlackList(gaIRanges_filtered, 
                blackList_list, verbose)
    }
    
    if(reportPDF){
        
        names_vec <- .computeExpNames(endoExoBamlist)
        exp_number <- length(names_vec)
        cols <- c("darkblue", "darkgreen", "darkred", "darkmagenta",
                "darkgray","darkorange", "darkcyan", "black",
                rainbow(exp_number-8))
        ssd_vec <- .computeSSD(ga, names_vec, verbose)
        gini_list <- .computeGini(gaIRanges, nbCPU, names_vec, verbose, 
                num_sim)
        estimated_thresholds <- .retrieveFilteringInfo(result_filtering, 
                .thresholdFUN)
        nb_total_reads <- .retrieveFilteringInfo(gaIRanges, .totalReadsFUN)
        nb_remaining_reads <- .retrieveFilteringInfo(gaIRanges_filtered, 
                .totalReadsFUN)
        perc_remaining_reads <- round((nb_remaining_reads*100)/nb_total_reads, 
                digits = 2)
        
        if(do_blacklist){
            nb_remaining_reads_blacklist <- 
                    .retrieveFilteringInfo(gaIRanges_filtered_blackList, 
                            .totalReadsFUN)
            nb_reads_inBlacklist <- nb_remaining_reads - 
                    nb_remaining_reads_blacklist
            perc_remaining_reads_blacklist <- 
                    round((nb_remaining_reads_blacklist*100)/nb_total_reads, 
                            digits = 2)
        }
        
        ## Building summary table
        table_summary <- cbind(gini_list[[1]], ssd_vec, estimated_thresholds,
                nb_total_reads, nb_remaining_reads, perc_remaining_reads)
        names_table <- c("Gini", "Gini.adjust", "SSD", "DupThres.",
                "TotalReads", "FilteredReads", "Perc.Kept.Reads")
        
        if(do_blacklist){
            table_summary <- cbind(table_summary, nb_reads_inBlacklist, 
                    nb_remaining_reads_blacklist, 
                    perc_remaining_reads_blacklist)
            names_table <- c(names_table, "ReadsInBlackList", 
                    "NbMinusBlackList", "PercMinusBlackList")
        }
        
        colnames(table_summary) <- names_table
        
        ## Evaluating the number of remaining reads according to increasing 
        ## duplicate threshold
        if(is.null(dupThreshold)){
            if(verbose)
                message("Testing the effect of different filtering thresholds")
            thresholdTestList <- .testTresholdEffect(gaIRanges, 
                    estimated_thresholds, verbose)
        }
        
        if(verbose)
            message("Building the PDF report")
        
        pdf(file=paste0(outputFolder, "duplicatesReport.pdf"), onefile = TRUE,
                ...)
        ## Page 1 - table Summary
        theme_p1 <- ttheme_minimal(base_size = 9)
        grid.table(table_summary, theme = theme_p1)
        
        ## Page 2 - Gini and Gini adj 
        colsp2 <- cols[seq_len(exp_number)]
        barplot(gini_list[[1]], beside = TRUE, col = colsp2)
        
        ## Page 3 - ssd
        barplot(ssd_vec, beside = TRUE, col = cols, xlab = "SSD", 
                axisnames = FALSE)
        
        ## Page 4 - plotting legend of p2 and p3
        plot.new()
        legend("topleft", names_vec, fill = colsp2)
        
        ## Page 5 - number of reads
        matReads <- rbind(nb_total_reads, nb_remaining_reads, 
                nb_remaining_reads_blacklist)
        colsp3 <- cols[seq_len(3)]
        barplot(matReads, beside = TRUE, col = colsp3, axisnames = FALSE,
                las = 2)
        text(seq(3, length(names_vec)*4, by = 4),  par("usr")[3]-0.25, adj = 1,
                xpd = TRUE, srt = 60, labels = names_vec)
        
        ## Page 6 - Plotting legend of p5
        categories_names <- c("Total", "Remaining")
        if(do_blacklist)
            categories_names <- c(categories_names, "MinusBlackList")
        plot.new()
        legend("topleft", categories_names, fill = colsp3)
        dev.off()
    }
    
    
    
    
    return(object)
}

# Note 1
# Each element of files_list contains two bam files:
# 1) the endogenous bam and 2) the exogenous bam
# The number of elements of the files_list is defined
# by the number of aligners that were used

setMethod(
        
        f = "filterDuplicates",
        
        signature = "ChIPSeqSpikeProcessList",
        
        definition = function(theObject, aligners = c("Bowtie", "Rsubread"), 
                nbCPU = 1, reportPDF = TRUE, outputFolder = NULL, num_sim = 1,
                dupThreshold = NULL, blackListEndo = NULL, blackListExo = NULL,
                verbose = TRUE, ...){
            
            datasetList(theObject) <- lapply(getDatasetList(theObject), 
                    function(object){
                        
                        bamPrep <- .prepareBAM(object, verbose)
                        groups <- rep(seq_len((length(bamPrep)-1)/2))
                        # see Note 1 above
                        files_list <- split(head(bamPrep, -1), groups)
                        
                        ## Filtering PCR duplicates and black lists
                        object <- .filterPCRduplicates(object, 
                                files_list, bamPrep[["paired"]], nbCPU, 
                                reportPDF, outputFolder, num_sim, dupThreshold,
                                blackListEndo, blackListExo, verbose)
                        
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


object1 <- getDatasetList(csdsProcessSE)[[1]]
endoBamBowtie(object1) <- ""
exoBamBowtie(object1) <- ""
endoBamRsubread(object1) <- ""
exoBamRsubread(object1) <- ""
object2 <- getDatasetList(csdsProcessSE)[[2]]
endoBamBowtie(object2) <- ""
exoBamBowtie(object2) <- ""
datasetList(csdsProcessSE) <- list(object1, object2)

filterDuplicates(csdsProcessSE)


