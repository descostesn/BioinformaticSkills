########################
# This file contains methods that were extracted (with modifications) from the
# htSeqTools package. Since this package will not be supported anymore from
# Bioconductor 3.8, the code below aims at making ChIPSeqSpike independent 
# from htSeqTools.
########################


#########################
# GENERICS DECLARATION
#########################

setGeneric(
        
        name = "ssdCoverage", 
        
        def = function(theObject){
            
            standardGeneric("ssdCoverage")
        },
        signature = "theObject")

setGeneric(
        
        name = "countRepeats", 
        
        def = function(theObject, counts){
            
            standardGeneric("countRepeats")
        },
        signature = "theObject")

setGeneric(
        
        name = "tabDuplReads", 
        
        def = function(theObject, nbCPU = 1){
            
            standardGeneric("tabDuplReads")
        },
        signature = "theObject")

setGeneric(
        
        name = "giniCoverage", 
        
        def = function(theObject, mc.cores = 1, numSim = 1){
            
            standardGeneric("giniCoverage")
        },
        signature = "theObject")

setGeneric(
        
        name = "filterDuplReads", 
        
        def = function(theObject, countrepeats, maxRepeats, mc.cores,
                fdrOverAmp = 0.01, negBinomUse = 0.999, components = 0,
                testThresholdVec = NULL, verbose = TRUE){
            
            standardGeneric("filterDuplReads")
        },
        signature = "theObject")


#########################
# METHODS
#########################

######################### common methods

setMethod(
        
        f = "countRepeats", 
        
        signature = "IRanges", 
        
        definition = function(theObject, counts) {
            
            tdf = data.frame(oldOrder=1:length(theObject),pos=start(theObject),
                    width=width(theObject), strand="+")
            tdf = tdf[order(tdf$pos, tdf$width),]
            if(length(unique(tdf$width)) == 1)
            {
                reps = Rle(tdf$pos)@lengths
            } else {
                reps = Rle(paste(tdf$pos, tdf$width, sep="."))@lengths
            }
            
            if(missing(counts)){
                return(reps)
            }else{
                cc <- as.list(reps)
                cc[reps>1] <- lapply(reps[reps>1], 
                        function(x) as.numeric(rep(x,x)))  
                readReps <- unlist(cc)[order(tdf$oldOrder)]
                ans <- list(reps=reps,readReps=readReps)
                return(ans)
            }
        }
)

######################### ssdCoverage

setMethod(
        
        f = "ssdCoverage", 
        
        signature = "GRanges", 
        
        definition = function(theObject) {
            
            covx <- coverage(theObject)
            cvs6 <- weighted.mean(sd(covx),w=sapply(covx, length))
            cvs6 <- 1000*cvs6/sqrt(length(theObject))
            return(cvs6)
        }
)

######################### tabDuplReads

setMethod(
        
        f = "tabDuplReads",
        
        signature = "IRangesList",
        
        definition = function(theObject, nbCPU = 1) {
            
            result <- mclapply(theObject, countRepeats, mc.cores = nbCPU)
            return(table(unlist(result)))
        }
)

######################### giniCoverage


.sampleRange<-function(x, chrLen, mc.cores){
    
    totReads<- sum(vapply(x,length,numeric(1)))
    chrReads<- totReads*(chrLen/sum(as.numeric(chrLen)))
    
    rangesl <- mcmapply(function(current_chrL, current_chrReads, 
                    current_IRanges){
                
                reads <- sample.int(current_chrL, as.integer(current_chrReads),
                        replace=T)
                minNb <- min(length(current_IRanges), 10000, current_chrL)
                len <- floor(mean(width(current_IRanges[1:minNb,])))
                ranges<-IRanges(start = reads, width = rep(len, length(reads)))
                return(ranges)
            }, chrLen, chrReads, x, mc.cores=mc.cores)
    
    return(rangesl)
}

.generateCounts<- function(sample, mc.cores){
    
    cove<- mclapply(sample, coverage, 
            mc.cores=mc.cores, mc.preschedule=FALSE)
    counts <- mclapply(cove, function(x){
                x <- Rle(values = runValue(x)[order(runValue(x))], 
                        lengths = runLength(x)[order(runValue(x))])
                result <- structure(array(runLength(x), dim=nrun(x), 
                                dimnames=structure(list(
                                                as.character(runValue(x))), 
                                        names="")), class="table")
            }, mc.cores=mc.cores, mc.preschedule=FALSE)
    return(counts)
}

.mergeChromoTable <- function(x){
    
    xnames <- unique(unlist(lapply(x, names)))
    xtable <- vector(length=length(xnames), mode="numeric")
    names(xtable) <- xnames
    
    for(i in xnames){
        
        for(j in names(x)){
            
            if(i %in% names(x[[j]])){
                xtable[i]<-xtable[i]+x[[j]][i]
            }
        }
    }
    return(xtable)
}

.plotLc <- function (x, general = FALSE, lwd = 2, xlab = "p", 
        ylab = "L(p)", main = "Lorenz curve", las = 1, ...){
    if (!general) 
        L <- x$L
    else L <- x$L.general
    plot(x$p, L, type = "l", main = main, lwd = lwd, xlab = xlab, 
            ylab = ylab, xaxs = "i", yaxs = "i", las = las, ...)
    abline(0, max(L))
}

.lorenzC <- function (x, n = rep(1, length(x)), plot = FALSE) 
{
    k <- length(x)
    o <- order(x)
    x <- x[o]
    n <- n[o]
    x <- n * x
    p <- c(0, cumsum(n)/sum(n))
    L <- c(0, cumsum(x)/sum(x))
    L2 <- L * mean(x)
    Lc <- list(p, L, L2)
    names(Lc) <- c("p", "L", "L.general")
    class(Lc) <- "Lc"
    
    if (plot)
        function() .plotLc(Lc)
    else
        return(Lc)
}

.ginifun <- function(indexes, nvec) {
    
    lorenz <- .lorenzC(indexes, nvec)
    gini<-1-sum(diff(lorenz$p) * (lorenz$L[-1] + lorenz$L[-length(lorenz$L)]))
    return(gini)
}

.plotProportion <- function(numvec, main){
    plot(log(numvec/sum(numvec)), type="h", main=main, 
            ylab='Proportion of bases (log)')
}

.plotRes <- function(x){
    
    main <- paste("Gini index: ", round(x[['gini.adjust']],4), sep="")
    
    #This line is recording the proportion plot. This can be displayed
    # by calling 'proportionFunc()'. The following line records the lorenz's
    proportionFUN <- function(){.plotProportion(x[["counts"]], main)}
    lorenzFUN <- .lorenzC(as.numeric(names(x[['counts']])), x[['counts']], 
            plot = TRUE)
    return(list(proportionFUN, lorenzFUN))
}

setMethod(
        
        f = "giniCoverage", 
        
        signature = "IRangesList",
        
        definition = function(theObject, mc.cores = 1, numSim = 1) {
            
            chrLengths<-unlist(lapply(theObject, 
                            function(y) max(end(ranges(y)))))
            
            # Simulating uniformily distributed data
            if(numSim>1) {
                message("\t\t Simulating uniformily distributed data ", numSim,
                        " times")
                simRange <- lapply(1:numSim, 
                        function(x) .sampleRange(theObject, chrLengths,
                                    mc.cores))
                simCounts <- lapply(simRange, .generateCounts, mc.cores)
                simCounts <- lapply(simCounts, .mergeChromoTable)
                simGini <- vapply(simCounts, 
                        function(x) .ginifun(as.numeric(names(x)), x), 
                        numeric(1))
                simGini <- unlist(simGini)
                simGinim <- mean(simGini)
                message("\t\t\t Average simulated gini: ", simGinim)
                message("\t\t\t Standard deviation: ", sd(simGini), "\n") 
                simGini <- simGinim
            } else {
                message("\t\t Simulating uniformily distributed data")
                simRange <- .sampleRange(theObject, chrLengths, mc.cores)
                simCounts <- .generateCounts(simRange, mc.cores)
                simCounts <- .mergeChromoTable(simCounts)
                simGini <- .ginifun(as.numeric(names(simCounts)), simCounts)
            }
            
            message("\t\t Calculating gini index of original data\n")
            counts <- .generateCounts(theObject, mc.cores)
            counts <- .mergeChromoTable(counts)
            gini <- .ginifun(as.numeric(names(counts)), counts)
            giniIdx <- list(gini= round(gini, digits = 2), 
                    gini.adjust= round(gini-simGini, digits = 2), 
                    counts = counts)
            
            ## Record the two functions (with specific values) to plot Lorenz
            ## and proportions later on
            plotFUN <- .plotRes(giniIdx)
            giVec <- unlist(giniIdx[c('gini', 'gini.adjust')])
            return(list(giniplots = plotFUN, ginicoeff = giVec))
        }
)




######################### filterDuplReads

.determineMaxRepeats <- function(countrepeats, negBinomUse, mc.cores, 
        fdrOverAmp){
    
    counts <- countrepeats
    use <- 1:sum(cumsum(counts)/sum(counts) < negBinomUse)
    if (length(use) <= 6) {
        components <- 1
    } else if (length(use) <= 8 & components > 2) {
        components <- 2
    } else if (length(use) <= 12 & components == 4) {
        components <- 3
    }
    fdr <- fdrEnrichedCounts(counts, use = use, 
            components = components, 
            mc.cores = mc.cores)$fdrEnriched
    maxRepeats <- max(5, match(TRUE, fdr < fdrOverAmp)-1)
    if (is.na(maxRepeats)) 
        maxRepeats <- max(as.numeric(counts)) + 1
    return(maxRepeats)
}

.filterObjectByMaxRepeats <- function(theObject, seqid, maxRepeats){
    return(mapply(function(x, id) x[id <= maxRepeats, ], theObject, seqid))
}

setMethod(
        
        f = "filterDuplReads",
        
        signature = "IRangesList",
        
        definition = function(theObject, countrepeats, maxRepeats, mc.cores,
                fdrOverAmp = 0.01, negBinomUse = 0.999, components = 0,
                testThresholdVec = NULL, verbose = TRUE) {
            
            if(verbose)
                message("\t\t Computing repeats")
            
            nrepeats <- lapply(theObject, countRepeats, countrepeats)
            seqid <- lapply(nrepeats, function(x) x[['readReps']])
            
            if(verbose)
                message("\t\t filtering")
            
            if(is.null(testThresholdVec)){
                
                if (is.null(maxRepeats)){
                    maxRepeats <- .determineMaxRepeats(countrepeats, 
                            negBinomUse, mc.cores, fdrOverAmp)
                }
                theObject <- .filterObjectByMaxRepeats(theObject, seqid, 
                        maxRepeats)
                
                return(list(filteredobject = theObject, 
                                threshold = maxRepeats))
            }else{
                
                filteredObjectList <- lapply(testThresholdVec, function(thres){
                            if(verbose)
                                message("\t\t\t Applying threshold ", thres)
                            return(.filterObjectByMaxRepeats(theObject, 
                                            seqid, thres))
                            
                        })
                return(filteredObjectList)
            }
            
        }
)


