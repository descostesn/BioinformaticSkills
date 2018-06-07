#############
# This script performs different statistics (parametric/non parametric) on a table contaning mean expression values that were used to construct boxplots.
# USE THE SCRIPT separate_files_fromMatrix.R BEFORE USING THIS ONE.
# Descostes April 2017
#############

library("NGSprofiling");
library("Rargs");
library("nortest")


##### !!!!!!!!!! code the parametric test if having normal distribution


################
# PARAMETERS
################

#parameters defined from the command line using RIO
paramsDefinition <- list();


paramsDefinition[["--valuesFileVec"]] <- list(variableName="values_file_vec", numeric=F, mandatory=T, description="Space separated list of files output by the script doing the boxplot considered or single matrix file with txt extension.");
paramsDefinition[["--outputFolder"]] <- list(variableName="output_folder", numeric=F, mandatory=T, description="Single character string defining the output folder.");
paramsDefinition[["--removeOutliers"]] <- list(variableName="remove_outliers", numeric=F, mandatory=T, description="Boolean indicating if outliers should be removed.", postConversion=as.logical);

paramsDefinition[["--expnamesVec"]] <- list(variableName="expnames_vec", numeric=F, mandatory=F, description="Space separated list of exp names corresponding to the exp defines in the previous parameter.", default = NA);

#values_file_vec <- c("/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/matrix-allpeaks_K27me3_allconditions-K27me3.WT.txt",
#                "/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/matrix-allpeaks_K27me3_allconditions-K27me3.Y358A.txt",
#                "/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/matrix-allpeaks_K27me3_allconditions-K27me3.F97A.txt",
#                "/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/matrix-allpeaks_K27me3_allconditions-K27me3.Y365A.txt",
#                "/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/matrix-allpeaks_K27me3_allconditions-K27me3.Trunc.txt")
#expnames_vec <- c("K27me3-WT", "K27me3-Y358A", "K27me3-F97A", "K27me3-Y365A", "K27me3-Trunc")
#output_folder <- c("/ifs/home/descon01/analysis/EED-mutant/boxplots/results/allpeaks/matrix_files/K27me3/with_outliers/")
#remove_outliers <- "FALSE";

#values_file_vec <- "/home/descostes/Documents/analysis/fact_ledgf/boxplots/rnaseqmay2018_comparison/based_upregulatedgenes_from_Day0ToDay6_downdKODay6/edgeR/TPM-matrix.txt"

################




##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);


if(!isTRUE(all.equal(length(output_folder), 1)))
{
	stop("output folder should be unique\n");
}


if(!isTRUE(all.equal(length(values_file_vec), length(expnames_vec))))
{
	stop("One exp name should be given per experiment\n");
}

checkingOutputFolder(output_folder);

#Reading the input table

extension <- unique(unlist(lapply(strsplit(basename(values_file_vec),"\\."),"[",2)));

if(isTRUE(all.equal(length(values_file_vec), 1) && isTRUE(all.equal(extension, "bed"))))
    stop("This script takes only one input matrix but not in bed format")


if(isTRUE(all.equal(extension, "bed")))
{
	list_exp <- lapply(values_file_vec, function(x){return(read.table(x, stringsAsFactors=F, header=TRUE))});
	list_exp <- lapply(list_exp, function(x){return(x$ratio)});
	
}else{
    if(isTRUE(all.equal(length(values_file_vec), 1))){
        
        current_mat <- read.table(values_file_vec, header=TRUE, stringsAsFactors = FALSE)
        expnames_vec <- colnames(current_mat)
        list_exp <- lapply(seq_len(ncol(current_mat)), function(i) current_mat[,i])
    }else{
        list_exp <- lapply(values_file_vec, function(x){return(read.table(x, stringsAsFactors=F, header=TRUE)[[1]])});
    }
}


#If defined, remove outliers

if(remove_outliers)
{
	list_exp <- lapply(list_exp, function(x){result_boxplot <- boxplot(x,plot=FALSE);
				
				if(!isTRUE(all.equal(length(result_boxplot$out),0))) return(x[-match(result_boxplot$out, x)])
				else return(x);})
}


####
# Part 1: Plotting the distribution of each list of values 
####

index_exp <- 1;
output_folder_distribution <- paste0(output_folder, "distribution/");
checkingOutputFolder(output_folder_distribution);

result <- lapply(list_exp, function(val_vec){
			
			# Plotting density without transformation
			png(filename=paste(output_folder_distribution, expnames_vec[index_exp], "distribution.png",sep=""), width = 600, height = 600, bg = "transparent");
			plot(density(val_vec));
			dev.off();
			
			#Plotting density with log transformation
			png(filename=paste(output_folder_distribution, expnames_vec[index_exp], "distribution-asinh.png",sep=""), width = 600, height = 600, bg = "transparent");
			plot(density(asinh(val_vec)), ylab="asinh(val)");
			dev.off();
			
			index_exp <<- index_exp + 1;
		});



####
# Part 2: Performing a mann whitney wilcoxon 2 by 2. 
####


output_folder_wilcox <- paste0(output_folder, "wilcoxon/");
checkingOutputFolder(output_folder_wilcox);

cat("Performing the mann whitney wilcoxon test on all combinations\n");

all_combination_vec <- combn(1:length(list_exp), 2);


result_wilcox_list <- apply(all_combination_vec, MARGIN=2, function(x){
			
			result <- wilcox.test(list_exp[[x[1]]], list_exp[[x[2]]], alternative="two.sided", mu=0, paired=FALSE, conf.int = TRUE);
			return(result);
		});


cat("output the results of the test\n");

names_test_vec <- apply(all_combination_vec, MARGIN=2, function(x){paste(expnames_vec[x],collapse="_")});
names(result_wilcox_list) <- names_test_vec;
pvalue_vec <- unlist(lapply(result_wilcox_list, "[[", "p.value"))
pvalue_bonferroni_vec <- p.adjust(pvalue_vec, method = "bonferroni")
pvalue_BH_vec <- p.adjust(pvalue_vec, method = "BH")

exact_pvalue_formated <- paste(names(pvalue_vec), ": ", pvalue_vec,  sep="")
bonferonni_pvalue_formated <- paste(names(pvalue_bonferroni_vec), ": ", pvalue_bonferroni_vec,  sep="")
BH_pvalue_formated <- paste(names(pvalue_BH_vec), ": ", pvalue_BH_vec,  sep="")

sink(paste0(output_folder_wilcox, "result_test-1.txt"));
print(result_wilcox_list);
print("exact p-val")
print(exact_pvalue_formated)
print("\n\n")
print("Bonferonni corrected p-val")
print(bonferonni_pvalue_formated)
print("\n\n")
print("BH corrected p-val")
print(BH_pvalue_formated)
sink();



####
# Part 3: Performing a Kruskal-Wallis test on the whole population (to avoid possible biases of the 2 by 2 comparison). Similar to a one-way anova for gaussian distribution 
####


output_folder_kruskal <- paste0(output_folder, "kruskal/");
checkingOutputFolder(output_folder_kruskal);

result_kruskal <- kruskal.test(list_exp);


sink(paste0(output_folder_kruskal, "result_test.txt"));
print(result_kruskal);
sink();




####
# Part 4: Verifying normality and performing a t test
####


cat("Computing t-test after verifying if the data are normally distributed\n");

output_folder_ttest <- paste0(output_folder, "ttest/");
checkingOutputFolder(output_folder_ttest);

# Performing a Kolmogorov-Smirnov test on each distribution

output_folder_AD <- paste0(output_folder_ttest, "Anderson-Darling/");
checkingOutputFolder(output_folder_AD);
index_exp <- 1;
result_AD_list <- lapply(list_exp, function(val_vec){
			
			result <- ad.test(val_vec);
			png(filename=paste(output_folder_AD, expnames_vec[index_exp],".png",sep=""), width = 600, height = 600, bg = "transparent")
			qqnorm(val_vec, main= paste0("AD test: ", result$p.val));
			qqline(val_vec);
			dev.off();
			index_exp <<- index_exp+1;
			return(NULL);
		});

# Performing the t test

result_ttest_list <- apply(all_combination_vec, MARGIN=2, function(x){
			
			#Testing if the variances are equal with a F-test
			result_ftest <- var.test(list_exp[[x[[1]]]], list_exp[[x[[2]]]], alternative="two.sided");
			
			result <- t.test(list_exp[[x[1]]], list_exp[[x[2]]], alternative="two.sided", mu= (mean(list_exp[[x[2]]]) - mean(list_exp[[x[1]]])), paired=FALSE, 
					var.equal = if(result_ftest$p.value < 0.05) FALSE else TRUE, 
					conf.level = 0.95);
			return(result);
		});

names(result_ttest_list) <- names_test_vec;

sink(paste0(output_folder_ttest, "result_test.txt"));
print(result_ttest_list);
sink();




