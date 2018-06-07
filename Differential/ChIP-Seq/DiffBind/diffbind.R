#################
## This script uses the bioconductor package DiffBind to perform differential 
## binding analysis with ChIP-Seq data. This script is written specifically and
## limited to the study of Gary, i.e, it handles only two conditions and the 
## contrasts are built on these two.
## -> Modifications enables to consider scaling factors in the context of 
## spiked-in data.
## Descostes May 2018
##
## SENS OF DB IS INVERTED ---> CORRECT
##
#################


library("DiffBind")
library("pryr")
library("Rargs")

!!!!!!!!!!!!!!! CORRECT THE INVERSION

##########################
## PARAMETERS
##########################
#
##cmd <- "/home/descostes/git/reinbergblab_sandbox/tools/differential_binding/diffbind/diffbind.R"


#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--csvFile"]] <- list(variableName="csv_file", numeric=F, mandatory=T, description="csv file containing the input parameter values")
paramsDefinition[["--analysisName"]] <- list(variableName="analysis_name", numeric=F, mandatory=T, description="String for the name of the analysis")
paramsDefinition[["--elongationSize"]] <- list(variableName="elongation_size", numeric=T, mandatory=T, description="Vector giving the elongation size for each exp. Should be of the form chip1, chip2, input1, input2")
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Path to the result folder (with slash)")
paramsDefinition[["--scalingFactors"]] <- list(variableName="scaling_factors", numeric=T, mandatory=T, description="Numeric vector giving the spiked-in factors. Should be NULL if none")

#Optional
paramsDefinition[["--minOverlap"]]=list(variableName="min_overlap", numeric= T, mandatory=F, description="Minimum nb of bp to overlapping to build the input matrix. Default: 0 (include all peaks)", default=0)
paramsDefinition[["--runParallel"]]=list(variableName="run_parallel", numeric=F, mandatory=F, description="If using all cpu available. Default = TRUE", postConversion=as.logical, default=TRUE)
paramsDefinition[["--analysisMethodVec"]]=list(variableName="analysis_method_vec", numeric=F, mandatory=F, description="This script only support DESEQ2 because using only one exp per condition", default="DESEQ2")
paramsDefinition[["--pValue"]]=list(variableName="p_value", numeric=T, mandatory=F, description="P-value to consider DB as significant. Default=0.05", default=0.05)
paramsDefinition[["--minQualThreshold"]]=list(variableName="min_qual_threshold", numeric=T, mandatory=F, description="Minimum sequencing quality for filtering reads. Default = 15", default=15)
paramsDefinition[["--verbose"]]=list(variableName="verbose", numeric=F, mandatory=F, description="Logical for printing processing messages. Default = TRUE", postConversion=as.logical, default=TRUE)
paramsDefinition[["--bpAroundSummit"]]=list(variableName="bp_around_summit", numeric=F, mandatory=F, description="FALSE or a positive value indicating the number of base pair around summits. Default: FALSE", postConversion=as.logical, default=FALSE)
paramsDefinition[["--nbReadsMinThres"]]=list(variableName="nb_reads_minThres", numeric=T, mandatory=F, description="If filter is a vector of values, ‘dba.count’ will return a vector of the same length, indicating how many intervals will be retained for each specified ‘filter’ level. Default: 0", default=0)
paramsDefinition[["--removeDuplicates"]]=list(variableName="remove_duplicates", numeric=F, mandatory=F, description="Logical indicating is duplicated should be removed. Default: TRUE", postConversion=as.logical, default=TRUE)
paramsDefinition[["--scaleControl"]]=list(variableName="scale_control", numeric=F, mandatory=F, description="See dba.count doc, changed default to FALSE", postConversion=as.logical, default=FALSE)
paramsDefinition[["--contrastColumns"]]=list(variableName="contrast_columns", numeric=F, mandatory=F, description="possible values (and can be combined) are: 'DBA_ID', 'DBA_TISSUE', 'DBA_FACTOR', 'DBA_CONDITION', 'DBA_TREATMENT', 'DBA_REPLICATE', 'DBA_CALLER'. Default: 'DBA_CONDITION'", default="DBA_CONDITION")

#csv_file <- "/ifs/home/descon01/analysis/fact_ledgf/diffbind/files_csv/ongenes/esc-polII.csv"
#analysis_name <- "esc-polII"
#elongation_size <- c(181,156,111,111)
#output_folder <- "/ifs/home/descon01/analysis/fact_ledgf/diffbind/results/ongenes/esc-polII/"
#scaling_factors <- c(0.0389700754636026, 0.0377913381497134)
#min_overlap <- 0
#run_parallel <- TRUE
#analysis_method_vec <- c("DESEQ2")  ## possible values are c("DESEQ2", "EDGER")
#p_value <- 0.05
#min_qual_threshold <- 15
#verbose <- TRUE
#annotation_GFF <- NULL   ## Annotations for counting reads on. If NULL, using the peak sets
#score_norm_vec <- c(DBA_SCORE_READS, DBA_SCORE_READS_FOLD, 
#        DBA_SCORE_READS_MINUS, DBA_SCORE_RPKM, DBA_SCORE_RPKM_FOLD, 
#        DBA_SCORE_TMM_READS_FULL, DBA_SCORE_TMM_READS_EFFECTIVE, 
#        DBA_SCORE_TMM_MINUS_FULL, DBA_SCORE_TMM_MINUS_EFFECTIVE, 
#        DBA_SCORE_TMM_READS_FULL_CPM, DBA_SCORE_TMM_READS_EFFECTIVE_CPM, 
#        DBA_SCORE_TMM_MINUS_FULL_CPM, DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM, 
#        DBA_SCORE_SUMMIT, DBA_SCORE_SUMMIT_ADJ, DBA_SCORE_SUMMIT_POS)
#bp_around_summit <- FALSE  ## FALSE or a positive value indicating the number of base pair
#nb_reads_minThres <- 0 ## If filter is a vector of values, ‘dba.count’ will return a vector of the same length, indicating how many intervals will be retained for each specified ‘filter’ level.
#remove_duplicates <- TRUE
#scale_control <- TRUE 
#filtering_function <- max ## Default is ‘max’, indicating that at least one sample should have a score of at least ‘filter’; other useful values include ‘sum’ (indicating that all the scores added together should be at
##                          ## least ‘filter’) and ‘mean’ (setting a minimum mean normalized count level). Users can supply their own function as well.
#contrast_columns <- "DBA_CONDITION"   ## possible values (and can be combined) are: "DBA_ID", "DBA_TISSUE", "DBA_FACTOR", "DBA_CONDITION", "DBA_TREATMENT", "DBA_REPLICATE", "DBA_CALLER"


#csv_file <- "/home/descostes/Documents/analysis/fact_ledgf/diffbind/esc_polII.csv"
#min_overlap <- 2 ## Default = 2
#run_parallel <- TRUE
#analysis_name <- "PolII-D0-D5"
#analysis_method_vec <- c("DESEQ2")  ## possible values are c("DESEQ2", "EDGER")
#p_value <- 0.05
#min_qual_threshold <- 15
#elongation_size <- c(181, 156, 111, 111)  ## Should be of the form chip1, chip2, input1, input2
#verbose <- TRUE
#output_folder <- "/home/descostes/Documents/analysis/fact_ledgf/diffbind/results/ESC/PolII_D0_vs_D5/"
#annotation_GFF <- NULL   ## Annotations for counting reads on. If NULL, using the peak sets
#
#score_norm_vec <- c(DBA_SCORE_READS, DBA_SCORE_READS_FOLD, 
#        DBA_SCORE_READS_MINUS, DBA_SCORE_RPKM, DBA_SCORE_RPKM_FOLD, 
#        DBA_SCORE_TMM_READS_FULL, DBA_SCORE_TMM_READS_EFFECTIVE, 
#        DBA_SCORE_TMM_MINUS_FULL, DBA_SCORE_TMM_MINUS_EFFECTIVE, 
#        DBA_SCORE_TMM_READS_FULL_CPM, DBA_SCORE_TMM_READS_EFFECTIVE_CPM, 
#        DBA_SCORE_TMM_MINUS_FULL_CPM, DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM, 
#        DBA_SCORE_SUMMIT, DBA_SCORE_SUMMIT_ADJ, DBA_SCORE_SUMMIT_POS)
#
#bp_around_summit <- FALSE  ## FALSE or a positive value indicating the number of base pair
#nb_reads_minThres <- 0 ## If filter is a vector of values, ‘dba.count’ will return a vector of the same length, indicating how many intervals will be retained for each specified ‘filter’ level.
#remove_duplicates <- TRUE
#scale_control <- TRUE 
#filtering_function <- max ## Default is ‘max’, indicating that at least one sample should have a score of at least ‘filter’; other useful values include ‘sum’ (indicating that all the scores added together should be at
#                          ## least ‘filter’) and ‘mean’ (setting a minimum mean normalized count level). Users can supply their own function as well.
#contrast_columns <- "DBA_CONDITION"   ## possible values (and can be combined) are: "DBA_ID", "DBA_TISSUE", "DBA_FACTOR", "DBA_CONDITION", "DBA_TREATMENT", "DBA_REPLICATE", "DBA_CALLER"
#scaling_factors <- c(0.0389700754636026, 0.0377913381497134)

#DBA_SCORE_READS                    raw read count for interval using only reads from ChIP                                                    
#DBA_SCORE_READS_FOLD               raw read count for interval from ChIP divided by read count for interval from control                     
#DBA_SCORE_READS_MINUS              raw read count for interval from ChIP minus read count for interval from control                          
#DBA_SCORE_RPKM                     RPKM for interval using only reads from ChIP                                                              
#DBA_SCORE_RPKM_FOLD                RPKM for interval from ChIP divided by RPKM for interval from control                                     
#DBA_SCORE_TMM_READS_FULL           TMM normalized (using edgeR), using ChIP read counts and Full Library size                                
#DBA_SCORE_TMM_READS_EFFECTIVE      TMM normalized (using edgeR), using ChIP read counts and Effective Library size                           
#DBA_SCORE_TMM_MINUS_FULL           TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size      
#DBA_SCORE_TMM_MINUS_EFFECTIVE      TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size 
#DBA_SCORE_TMM_READS_FULL_CPM       same as ‘DBA_SCORE_TMM_READS_FULL’, but reported in counts-per-million.                                   
#DBA_SCORE_TMM_READS_EFFECTIVE_CPM  same as ‘DBA_SCORE_TMM_READS_EFFECTIVE’, but reported in counts-per-million.                              
#DBA_SCORE_TMM_MINUS_FULL_CPM       same as ‘DBA_SCORE_TMM_MINUS_FULL’, but reported in counts-per-million.                                   
#DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM  same as ‘DBA_SCORE_TMM_MINUS_EFFECTIVE’, but reported in counts-per-million.                              
#DBA_SCORE_SUMMIT                   summit height (maximum read pileup value)                                                                 
#DBA_SCORE_SUMMIT_ADJ               summit height (maximum read pileup value), normalized to relative library size                            
#DBA_SCORE_SUMMIT_POS               summit position (location of maximum read pileup)                               

#DO-BRD4_mm10: 0.0490711229003939
#D13-BRD4_mm10: 0.0352673389491398
#D5-RNAPolII_mm10: 0.0377913381497134
#D5-BRD4_mm10: 0.0424813040842758
#D0-RNAPolII_mm10: 0.0389700754636026
#D13-RNAPolII_mm10: 0.0344027544221301
#Day0_HDGF2_mm10: 0.033213175812983
#Day13_HDGF2_mm10: 0.0336824008950088
#Day0-input_mm10: 0.0280413933588895
#Day5-SPT16_mm10: 0.0322579771074169
#Day13-SPT16_mm10: 0.032329816894586
#Day5_input_biorep1_mm10: 0.0295097500952058
#Day0-SPT16_mm10: 0.0407366257380356
#Day5_HDGF2_mm10: 0.0302535959272851
#Day13_input_biorep1_mm10: 0.0258990089899604


#Warning messages:
#        1: In if (bParallel) { :
#                    the condition has length > 1 and only the first element will be used
#            2: Some groups have no replicates. Results may be unreliable. 
#            3: In if (DBA$config$parallelPackage == DBA_PARALLEL_MULTICORE) { :
#                        the condition has length > 1 and only the first element will be used
#                4: In checkForExperimentalReplicates(object, modelMatrix) : 
#                        
#                        Deprectation note: Analysis of designs without replicates will be removed
#                in the Oct 2018 release: DESeq2 v1.22.0, after which DESeq2 will give an error.
#                
#                5: In checkForExperimentalReplicates(object, modelMatrix) : 
#                        
#                        The design matrix has the same number of samples and coefficients to fit,
#                estimating dispersion by treating samples as replicates. This analysis
#                is not useful for accurate differential expression analysis, and arguably
#                not for data exploration either, as large differences appear as high dispersion.
                



#########################
# FUNCTION
#########################


.global2Str <- function(x){force(x); return(as.character(substitute(x)))}



continueAnalysis <- function(dba_struct, outputFold, contrastCols, condVec){
    
    if(!file.exists(outputFold))
        dir.create(outputFold, recursive = TRUE)
    
    pdf(file = paste0(outputFold, "/correlation_heatmap.pdf"))
    plot(dba_struct)
    dev.off()
    
    if(verbose)
        message("\t\t\t Establishing contrast on: ", 
                paste(contrastCols, collapse = " "))
    
    dba_struct <-  dba.contrast(DBA = dba_struct, 
            group1 = dba_struct$masks[[condVec[1]]],
            group2 = dba_struct$masks[[condVec[2]]],
            categories= sapply(contrastCols, get))
    
    if(verbose)
        message("\t\t\t #### Performing differential analysis")
    
    for(subInput in c(FALSE, TRUE)){
        
        for(fullLibrarySize in c(FALSE, TRUE)){
            
            if(verbose)
                message("\t\t\t\t ", if(subInput) "with " else "without ", 
                        "input subtraction and ", 
                        if(fullLibrarySize) "library size " 
            else "reads in peaks size ", "normalization")
            
            outputFold2 <- paste0(outputFold, 
                    if(subInput) "/inputSub" else "/noInputSub",
                    if(fullLibrarySize) "/librarySizeNorm" 
                            else "/readsInPeaksNorm")
            
            if(!file.exists(outputFold2))
                dir.create(outputFold2, recursive = TRUE)
            
            ## bTagwise is false because script considers only one rep
            dba_struct <- dba.analyze(DBA = dba_struct,
                    bSubControl = subInput, 
                    bFullLibrarySize = fullLibrarySize,
                    bTagwise = FALSE)
            
            if(verbose)
                message("\t\t\t\t\t Retrieving differential report")
            
            dba_struct.DB2FOLD <- dba.report(DBA = dba_struct, 
                    fold = 2, 
                    bCounts = TRUE,
                    file = "2fold",
                    initString = paste0(outputFold2, "/report"))
            
            dba_struct.DB2FOLDUP <- dba_struct.DB2FOLD[
                    dba_struct.DB2FOLD$Fold > 0, ]
            
            dba_struct.DB2FOLDDown <- dba_struct.DB2FOLD[
                    dba_struct.DB2FOLD$Fold < 0, ]
            
            dba_struct.ALLFOLD <- dba.report(DBA = dba_struct, 
                    fold = 0, 
                    bCounts = TRUE,
                    file = "all",
                    initString = paste0(outputFold2, "/report"))
            
            dba_struct.ALLFOLDUP <- dba_struct.ALLFOLD[
                    dba_struct.ALLFOLD$Fold > 0, ]
            
            dba_struct.ALLFOLDDown <- dba_struct.ALLFOLD[
                    dba_struct.ALLFOLD$Fold < 0, ]
            
            
            summary_str <- c(paste0("p-val < ", unique(dba_struct$config$th)),
                    paste0(" - Nb 2 fold: ", length(dba_struct.DB2FOLD)),
                    paste0(" - Nb 2 fold up: ", length(dba_struct.DB2FOLDUP)),
                    paste0(" - Nb 2 fold down: ", 
                            length(dba_struct.DB2FOLDDown)),
                    paste0(" - All DB: ", length(dba_struct.ALLFOLD)),
                    paste0(" - All DB up: ", length(dba_struct.ALLFOLDUP)),
                    paste0(" - All DB down: ", length(dba_struct.ALLFOLDDown)))
            
            dba_list <- list(dba_struct.ALLFOLD, dba_struct.ALLFOLDDown, 
                    dba_struct.ALLFOLDUP, dba_struct.DB2FOLD, 
                    dba_struct.DB2FOLDDown, dba_struct.DB2FOLDUP)
            gff_namesList <- list("ALLFOLD", "ALLFOLDDown", 
                    "ALLFOLDUP", "DB2FOLD", 
                    "DB2FOLDDown", "DB2FOLDUP")
            
            invisible(mapply(function(dba, name){
                                if(!is.null(dba) && 
                                        !isTRUE(all.equal(length(dba), 0)))
                                    rtracklayer::export(dba, 
                                            con=paste0(outputFold2, "/", 
                                                    name, ".gff"), 
                                            format="GFF")
                            }, dba_list, gff_namesList))
            
            
            write(summary_str, file = paste0(outputFold2, "/summary.txt"),
                    ncolumns = 1)
            
            enoughDB <- all(sapply(list(dba_struct.ALLFOLD, 
                                    dba_struct.ALLFOLDUP, 
                                    dba_struct.ALLFOLDDown), 
                            function(x) length(x) > 1))
            
            if(enoughDB){
                
                if(verbose)
                    message("\t\t\t\t\t Plotting correlation on DB peaks")
                pdf(file = paste0(outputFold2, "/DBheatmap.pdf"), width = 10, 
                        height = 10)
                plot(dba_struct, contrast = 1)
                dev.off()
                
#                if(verbose)
#                    message("\t\t\t\t\t Plotting volcano")
#                pdf(file = paste0(outputFold2, "/volcano.pdf"))
#                dba.plotVolcano(dba_struct)
#                dev.off()
                
#                if(verbose)
#                    message("\t\t\t\t\t Plotting cluster")
#                pdf(file = paste0(outputFold2, "/cluster.pdf"), width = 10, 
#                        height = 10)
#                dba.plotHeatmap(dba_struct, contrast=1, correlations = FALSE)
#                dev.off()
            }else{
                message("\t\t\t\t\t /-- No DB genes found --/")
            }
            
            
            if(verbose)
                message("\t\t\t\t\t Plotting PCA")
            pdf(file = paste0(outputFold2, "/PCA.pdf"))
            dba.plotPCA(dba_struct, DBA_TISSUE, label = DBA_CONDITION)
            dev.off()
            
            if(verbose)
                message("\t\t\t\t\t Plotting MA")
            pdf(file = paste0(outputFold2, "/MA.pdf"))
            dba.plotMA(dba_struct)
            abline(h = 2, col = "red")
            abline(h = -2, col = "red")
            dev.off()
            
        }
    }
    
    
}


#########################
# MAIN
#########################


# Retreives the parameters
getParams(paramsDefinition);



## Can be included later as param

filtering_function <- max ## Default is ‘max’, indicating that at least one sample should have a score of at least ‘filter’; other useful values include ‘sum’ (indicating that all the scores added together should be at
                          ## least ‘filter’) and ‘mean’ (setting a minimum mean normalized count level). Users can supply their own function as well.

score_norm_vec <- c(DBA_SCORE_READS, DBA_SCORE_READS_FOLD,
        DBA_SCORE_READS_MINUS, DBA_SCORE_RPKM, DBA_SCORE_RPKM_FOLD,
        DBA_SCORE_TMM_READS_FULL, DBA_SCORE_TMM_READS_EFFECTIVE,
        DBA_SCORE_TMM_MINUS_FULL, DBA_SCORE_TMM_MINUS_EFFECTIVE,
        DBA_SCORE_TMM_READS_FULL_CPM, DBA_SCORE_TMM_READS_EFFECTIVE_CPM,
        DBA_SCORE_TMM_MINUS_FULL_CPM, DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM,
        DBA_SCORE_SUMMIT, DBA_SCORE_SUMMIT_ADJ, DBA_SCORE_SUMMIT_POS)

annotation_GFF <- NULL

## Verify parameters

if(!all(analysis_method_vec %in% c("DESEQ2", "EDGER")))
    stop("Analysis methods should be DESEQ2 or EDGER")

if(!(isTRUE(all.equal(length(csv_file), 1)) || 
            isTRUE(all.equal(file_ext(csv_file), 1))))
    stop("The sample sheet should be unique and in csv format")

csv_validate <- read.csv(csv_file)

if(!(isTRUE(all.equal(length(unique(csv_validate$PeakCaller)), 1)) ||
            isTRUE(all.equal(unique(csv_validate$PeakCaller), "narrow"))))
    stop("Use the narrow value in PeakCaller column")

peakformat <- unique(csv_validate$PeakFormat)
if(!(isTRUE(all.equal(length(peakformat), 1)) ||
            (isTRUE(all.equal(peakformat, "broadPeak")) &&
                                isTRUE(all.equal(peakformat, "narrowPeak")))))
                    stop("Use the same peakFormat and use macs2 .narrowPeak",
                            " or .broadPeak files")

if(!file.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)

if(!(is.null(annotation_GFF) || (isTRUE(all.equal(length(annotation_GFF), 1)) 
                && isTRUE(all.equal(file_ext(annotation_GFF), "gff")))))
    stop("Annotation file should be unique and in .gff format")

## TO REMOVE AFTER COMPLETING CODE
if(!is.null(annotation_GFF))
    stop("The code is not adapted yet to handle a user defined annotation ",
            "file. The gff should be read in a GRange object given to the ",
            "same variable")

condition_vec <- sort(as.character(unique(csv_validate$Condition)))

if(!isTRUE(all.equal(length(condition_vec), 2)))
    stop("This script is designed to handle only two conditions")

if(!isTRUE(all.equal(condition_vec, c("condition1", "condition2"))))
    stop("For this script, allowed values for condition are 'condition1 ",
            "and condition2")

if(!isTRUE(all.equal(nrow(csv_validate), 2)))
    stop("This script is made to run without biological replicates")

if(isTRUE(all.equal(nrow(csv_validate), 2)) && 
        any(analysis_method_vec %in% "EDGER"))
    stop("DiffBind does not allow to use edgeR without replicates")


## PART 1: PERFORMING THE DIFFERENTIAL BINDING ANALYSIS

if(verbose){
    message(rep("#", 10))
    message("## PERFORMING THE DIFFERENTIAL BINDING ANALYSIS")
    message(rep("#", 10))
}


for(analysis_method in analysis_method_vec){
    
    if(verbose)
        message("# Method: ", analysis_method)
    
    
    amethod <- isTRUE(all.equal(analysis_method, "DESEQ2"))
    amethod <- if(amethod) DBA_DESEQ2 else DBA_EDGER
    
    if(verbose)
        message("\t ## Building dba object")
    
    ##Building initial object
    dba_obj <- dba(minOverlap = min_overlap, 
            sampleSheet = csv_file,
            config = data.frame(RunParallel = run_parallel,
                    reportInit = analysis_name,
                    DataType = DBA_DATA_GRANGES,
                    AnalysisMethod = amethod, 
                    minQCth = min_qual_threshold, 
                    fragmentSize = elongation_size, 
                    th = p_value, 
                    bUsePval = TRUE), scoreCol = 5)
    
    number_report <- data.frame(totalPeaks = nrow(dba_obj$merged),
            commonPeaks = nrow(dba_obj$binding),
            peaks1 = as.vector(dba.show(dba_obj)$Intervals[[1]]),
            peaks2 = as.vector(dba.show(dba_obj)$Intervals[[2]]))
    
    if(verbose){
        message("Initial figures: ", print(number_report))
        message("\t\t Generating correlation heatmap")
    }
    
    output_folder_tmp <- paste0(output_folder, analysis_method)
    if(!file.exists(output_folder_tmp))
        dir.create(output_folder_tmp, recursive = TRUE)
    pdf(file = paste0(output_folder_tmp, "/correlation_heatmap.pdf"))
    plot(dba_obj)
    dev.off()
    
    for(score_norm in score_norm_vec){
        
        score_method_name <- switch(as.character(score_norm),
                "1" =.global2Str(DBA_SCORE_RPKM),
                "2" = .global2Str(DBA_SCORE_RPKM_FOLD),
                "3" = .global2Str(DBA_SCORE_READS), 
                "4" = .global2Str(DBA_SCORE_READS_FOLD),
                "5" = .global2Str(DBA_SCORE_READS_MINUS),
                "6" = .global2Str(DBA_SCORE_TMM_MINUS_FULL),
                "7" = .global2Str(DBA_SCORE_TMM_MINUS_EFFECTIVE),
                "8" = .global2Str(DBA_SCORE_TMM_READS_FULL),
                "9" = .global2Str(DBA_SCORE_TMM_READS_EFFECTIVE), 
                "10" = .global2Str(DBA_SCORE_TMM_MINUS_FULL_CPM),
                "11" = .global2Str(DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM),
                "12" = .global2Str(DBA_SCORE_TMM_READS_FULL_CPM),
                "13" = .global2Str(DBA_SCORE_TMM_READS_EFFECTIVE_CPM),
                "101" = .global2Str(DBA_SCORE_SUMMIT),
                "102" = .global2Str(DBA_SCORE_SUMMIT_ADJ),
                "103" = .global2Str(DBA_SCORE_SUMMIT_POS))
        
        if(verbose)
            message("\t\t ### Using ", score_method_name)
        output_folder_tmp2 <- paste0(output_folder_tmp, "/", 
                score_method_name)
        
        
        if(verbose)
            message("\t\t\t Computing counts")
        
        if(!is.null(annotation_GFF)){
            
            dba_obj2 <- dba.count(DBA = dba_obj,
                    peaks = annotation_GFF, 
                    minOverlap = min_overlap, 
                    score = score_norm, 
                    bLog = if(isTRUE(all.equal(score_norm, 2)) || 
                                    isTRUE(all.equal(score_norm, 4))) 
                                TRUE else FALSE,
                    summits = bp_around_summit,
                    filter= nb_reads_minThres, 
                    bRemoveDuplicates = remove_duplicates,
                    bScaleControl = scale_control,
                    filterFun = filtering_function, 
                    readFormat = DBA_READS_BAM)
            }else{
                
                dba_obj2 <- dba.count(DBA = dba_obj,
                        minOverlap = min_overlap, 
                        score = score_norm, 
                        bLog = if(isTRUE(all.equal(score_norm, 2)) || 
                                        isTRUE(all.equal(score_norm, 4))) TRUE
                        else FALSE,
                        summits = bp_around_summit,
                        filter= nb_reads_minThres, 
                        bRemoveDuplicates = remove_duplicates,
                        bScaleControl = scale_control,
                        filterFun = filtering_function, 
                        readFormat = DBA_READS_BAM)
        }
        
        
        if(!is.null(scaling_factors) && 
                isTRUE(all.equal(score_method_name, "DBA_SCORE_READS"))){
            
            for(isspiked in c(TRUE, FALSE)){
                
                if(verbose)
                    message("\t\t\t ###--------- ", 
                            if(isspiked) "spiked" else "no spike")
                
                if(isspiked){ ## Applying scaling factors
                    
                    dba_obj_spiked <- dba_obj2
                    dba_obj_spiked$peaks <- mapply( 
                            function(peaks_df, sfact){
                                peaks_df$Score <- peaks_df$Score * sfact
                                return(peaks_df)
                            }, dba_obj_spiked$peaks, scaling_factors, 
                            SIMPLIFY = FALSE)
                    
                    dba_obj_spiked$binding <- cbind(dba_obj2$binding[,1:3], 
                            sweep(dba_obj2$binding[,4:5], 
                                    2, scaling_factors, "*"))
                    output_folder_tmp3 <- paste0(output_folder_tmp2, "/spiked")
                            
                    continueAnalysis(dba_obj_spiked, output_folder_tmp3, 
                            contrast_columns, condition_vec)
                    
                }else{ ## Continue without spike-in
                    continueAnalysis(dba_obj2, output_folder_tmp2, 
                            contrast_columns, condition_vec)
                }
            }
        }else{
            continueAnalysis(dba_obj2, output_folder_tmp2, contrast_columns, 
                    condition_vec)
        }
    }
}



## PART 2: BUILDING REPORT

if(verbose){
    message(rep("#", 10))
    message("## BUILDING REPORTS")
    message(rep("#", 10))
    message("Retrieving summary report files")
}

report_vec <- list.files(path = substr(output_folder, 1, 
                nchar(output_folder) - 1), pattern = "summary.txt", 
        full.names = TRUE, recursive = TRUE)
names_vec <- dirname(sub(paste0(output_folder, "DESEQ2/"), "", report_vec))


if(verbose)
    message("Reading reports into a list")
report_list <- lapply(report_vec, partial(readLines))
names(report_list) <- names_vec

if(verbose)
    message("Retrieving: Max DB, max DB 2 folds, max up, max down, max up ",
            "2 folds, max down 2 folds")

retrieveMaxIndex <- function(indexNum){
    
    result <- unlist(lapply(report_list, function(x){
                                as.numeric(strsplit(x[indexNum], ":")[[1]][2])
                            }))
    indexes <- which(result == max(result))
    return(as.numeric(indexes))
}


indexes_list <- lapply(c(5, 2, 6, 7, 3, 4), retrieveMaxIndex)
names(indexes_list) <- c("MaxDB", "MaxDB2folds", "MaxUp", "MaxDown", 
        "MaxUp2Folds", "MaxDown2Folds")
names_report_vec <- c("pvalThres", "2folds", "2foldsUp", "2foldsdown", "all", 
        "allup", "alldown")

matrices_list <- lapply(indexes_list, function(index_vec){
            
            mat <- do.call(rbind, report_list[index_vec])
            colnames(mat) <- names_report_vec
            
            if(!is.vector(mat[,-1]))
                mat <- apply(mat[,-1], 2, function(x) 
                        as.numeric(strsplit(x, ":")[[1]][2]))
            else{
                mat <- as.numeric(unlist(lapply(strsplit(mat[,-1], ":"),
                                        "[",2)))
                names(mat) <- names_report_vec[-1]
            }
            return(mat)
        })

mat_val <- do.call(rbind, matrices_list)

parameters_list <- lapply(indexes_list, function(index_vec) 
            names(report_list[index_vec]))
parameters_vec <- sapply(parameters_list, partial(paste), collapse="--")

write.table(mat_val, file=paste0(output_folder, "report.txt"), quote=FALSE, 
        sep="\t", col.names=TRUE, row.names=TRUE)
write.table(parameters_vec, file=paste0(output_folder, "report.txt"), 
        quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, append = TRUE)
write.table(number_report, file=paste0(output_folder, "report.txt"), 
        quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, append = TRUE)
