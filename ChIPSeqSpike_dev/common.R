.section <- function(msg, verbose, repNB){
    
    if(verbose){
        headstr <- "###################################"
        message(rep(headstr, repNB))
        message(msg)
        message(rep(headstr, repNB))
        message("\n")
    }
}

.call_fun <- function(f, ...) f(...)


.prepareBAM <- function(object, verbose){
    
    retrieveFastq <- list(getFASTQ1, getFASTQ2)
    fqVec <- lapply(retrieveFastq, .call_fun, object)
    paired <- !isTRUE(all.equal(length(fqVec[[2]]), 0))
    
    msg <- paste0("Processing ", fqVec[[1]], if(paired) paste0("and ", 
                        fqVec[[2]]), " alignment")
    .section(msg, verbose, 3)
    
    retrieveBam <- list(getEndoBAMRsubread, getExoBAMRsubread, 
            getEndoBAMBowtie, getExoBAMBowtie)
    result <- vector("list", 5)
    result <- lapply(retrieveBam, .call_fun, object)
    names_vec <- c(vapply(result, function(x) file_path_sans_ext(basename(x)),
                    character(1)), "paired")
    result <- c(result, paired)
    names(result) <- names_vec
    
    if(isTRUE(all.equal(length(result[[1]]), 0)) &&
            isTRUE(all.equal(length(result[[2]]), 0)))
        stop("No BAM aligned to ", getGenome(object),
                " can be found")
    return(result)
}


.endoNotExoError <- function(object, alignment_source){
    
    stop("For ", getFASTQ1(object), ", reads, ",
            "they were ", alignment_source, " aligned to ", 
            getGenome(object), " but not to ", 
            getExoGenome(object))
}

