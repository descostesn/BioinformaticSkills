###############
# This script performs ma plot from given objects generating all combinations and scaling according to the largest values on Y-axis.
# March 2018
###############

library("edgeR")
library("Rargs")



#########################
# PARAMETERS
#########################


#parameters defined from the command line using RIO
paramsDefinition <- list();

#firstly, definition of mandatory parameters that have to be given by the user in the command line

paramsDefinition[["--objectVec"]] <- list(variableName="object_vec", numeric=F, mandatory=T, description="Space separated list of objects of which mean values will be used");
paramsDefinition[["--expnameVec"]] <- list(variableName="expname_vec", numeric=F, mandatory=T, description="Corresponding experiment names");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Folder to which results will be written");
paramsDefinition[["--xAxisMultiplyingFactor"]] <- list(variableName="x_axis_multiplying_factor", numeric=T, mandatory=T, description="Scaling factor to change the x-axis");
paramsDefinition[["--xLabel"]] <- list(variableName="x_label", numeric=F, mandatory=T, description="Label to display on the x-axis");

paramsDefinition[["--xAxisTicks"]] <- list(variableName="x_axis_ticks", numeric=T, mandatory=F, description="Vector of numeric ticks for the x-axis (will also define the min and max.", default=NA);
paramsDefinition[["--yAxisTicks"]] <- list(variableName="y_axis_ticks", numeric=T, mandatory=F, description="Vector of numeric ticks for the y-axis (will also define the min and max.", default=NA);


#object_vec <- c("/ifs/home/descon01/analysis/EED-mutant/objects/for_revision/groupI_strongpeaks/spNA/GenesSelect_E14_WT_H3K27me3_spaceSizes_NA_ProLength_15000",
#                "/ifs/home/descon01/analysis/EED-mutant/objects/for_revision/groupI_strongpeaks/spNA/GenesSelect_E14_Y365A_H3K27me3_spaceSizes_NA_ProLength_15000",
#                "/ifs/home/descon01/analysis/EED-mutant/objects/for_revision/groupI_strongpeaks/spNA/GenesSelect_E14_F97A_H3K27me3_spaceSizes_NA_ProLength_15000",
#                "/ifs/home/descon01/analysis/EED-mutant/objects/for_revision/groupI_strongpeaks/spNA/GenesSelect_E14_Trunc_H3K27me3_spaceSizes_NA_ProLength_15000",
#                "/ifs/home/descon01/analysis/EED-mutant/objects/for_revision/groupI_strongpeaks/spNA/GenesSelect_E14_Y358A_H3K27me3_spaceSizes_NA_ProLength_15000")
#
#expname_vec <- c("WT_H3K27me3",
#        "Y365A_H3K27me3",
#        "F97A_H3K27me3",
#        "Trunc_H3K27me3", 
#        "Y358A_H3K27me3")
#
#
#output_folder <- c("/ifs/home/descon01/analysis/EED-mutant/maplot/for_revision/strongpeaks/")
#x_axis_multiplying_factor <- 10000
#x_label <- "Relative bin counts"
#
#x_axis_ticks <- c(5,10,15,20)
#y_axis_ticks <- c(-15, -10, -5, 0, 5)




#########################
# MAIN
#########################


# Retreives the parameters
getParams(paramsDefinition);


if(!file.exists(output_folder))
    dir.create(output_folder, recursive=TRUE)

if(!isTRUE(all.equal(length(object_vec), length(expname_vec))))
    stop("One expname should be given for each object")

if(!isTRUE(all.equal(length(output_folder),1)))
    stop("This script considers only one output folder")


message("Retrieving the objects")

object_list <- lapply(object_vec, function(x){
            message("Loading ", x)
            load(x)
            current <- genesSelected
            return(current)
        })


message("Retrieving count values inside")

count_inside_list <- lapply(object_list, function(x){
            values_inside <- lapply(x, "[[", "probe.valueINSIDE")
            sum_values <- unlist(lapply(values_inside, sum))
            return(sum_values)})
rm(object_list)
gc()
result_matrix <- do.call(cbind, count_inside_list)
rm(count_inside_list)
gc()


#Creating a list of matrices for each combination and filtering them properly

matrix_combn_index <- cbind(combn(seq_len(ncol(result_matrix)), 2), combn(rev(seq_len(ncol(result_matrix))), 2))
result_matrix_list <- apply(matrix_combn_index, 2, function(numcols, resultmat){
            
            tmp_mat <- cbind(resultmat[,numcols[1]], resultmat[,numcols[2]])
            
            message("Processing ", expname_vec[numcols[1]], "-", expname_vec[numcols[2]])
            ## Removing null and negative values
            exp_keep_indexes <- intersect(which(tmp_mat[,1] > 0), which(tmp_mat[,2] > 0))
            tmp_mat <- tmp_mat[exp_keep_indexes,]
            message("\t The number of kept values is: ", length(exp_keep_indexes))
            return(tmp_mat)
        }, result_matrix)


## Performing MA plot for all combinations

log2fc_list <- lapply(result_matrix_list, 
        function(current_mat){
            return(log2((current_mat[,2]/current_mat[,1])))})
name_output_list <- apply(matrix_combn_index, 2, function(x){return(paste(expname_vec[x[1]], expname_vec[x[2]], sep="-"))})

if(is.na(y_axis_ticks[1])){
    max_ylab <- max(unlist(lapply(log2fc_list,max)))
    min_ylab <- min(unlist(lapply(log2fc_list,min)))
}

if(is.na(x_axis_ticks[1])){
    min_xlab <- min(unlist(lapply(result_matrix_list, function(x){return(min(log2(x*x_axis_multiplying_factor)))})))
    max_xlab <- max(unlist(lapply(result_matrix_list, function(x){return(max(log2(x*x_axis_multiplying_factor)))})))
}


message("Generating MA plot")
mapply(function(mat, log2fc, nameoutput){
            pdf(file=paste0(output_folder, nameoutput, ".pdf"))
            maPlot(mat[,1]*x_axis_multiplying_factor, mat[,2]*x_axis_multiplying_factor, logFC=log2fc, normalize=FALSE, de.tags=NULL, smooth.scatter=TRUE, lowess=TRUE, ylab="log2(FC)", xlab= x_label, 
                    xaxt= if(!is.na(x_axis_ticks[1])) "n",
                    yaxt= if(!is.na(y_axis_ticks[1])) "n",
                    xlim= if(!is.na(x_axis_ticks[1])) c(min(x_axis_ticks), max(x_axis_ticks)) else c(min_xlab, max_xlab),
                    ylim= if(!is.na(y_axis_ticks[1])) c(min(y_axis_ticks), max(y_axis_ticks)) else c(min_ylab, max_ylab))
            if(!is.na(x_axis_ticks[1]))
                axis(1, xaxp=c(min(x_axis_ticks), max(x_axis_ticks), length(x_axis_ticks)-1), las=2)
            if(!is.na(y_axis_ticks[1]))
                axis(2, xaxp=c(min(y_axis_ticks), max(y_axis_ticks), length(y_axis_ticks)-1), las=2)
            dev.off()
        }, result_matrix_list, log2fc_list, name_output_list)

