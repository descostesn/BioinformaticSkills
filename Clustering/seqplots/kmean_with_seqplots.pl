#!/usr/bin/env perl
# This script performs heatmap of kmean clustering of chip-seq data using annotations.
# It uses the input file files_X.conf of format: BIGWIG;BIGWIGNAMEVEC;BEDVEC;BEDNAMEVEC;ORGANISM;GENOMEVERSION;BINSIZE;PROGILELENGTHBEFORE;PROFILELENGTHAFTER;TYPEVALUE;OUTPUTFOLDER;NBOFGROUPS;INCLUDEEXPVEC;clustering_method;sort_rows;autoscale;zmin;zmax  with optional color vec
#Nicolas Descostes feb 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 18 && scalar(@arguments_tab) != 19)
{
	die("Missing arguments for multiple_techniques.pl: BIGWIG;BIGWIGNAMEVEC;BEDVEC;BEDNAMEVEC;ORGANISM;GENOMEVERSION;BINSIZE;PROGILELENGTHBEFORE;PROFILELENGTHAFTER;TYPEVALUE;OUTPUTFOLDER;NBOFGROUPS;INCLUDEEXPVEC;clustering_method;sort_rows;autoscale;zmin;zmax with optional color vec\n\n");
}


my $bigwig_vec = $arguments_tab[0];
my $bigwig_name_vec = $arguments_tab[1];
my $bed_vec = $arguments_tab[2];
my $bed_name_vec = $arguments_tab[3];
my $organism = $arguments_tab[4];
my $genome_version = $arguments_tab[5];
my $bin_size = $arguments_tab[6];
my $profile_length_before = $arguments_tab[7];
my $profile_length_after = $arguments_tab[8];
my $type_value = $arguments_tab[9];
my $output_folder = $arguments_tab[10];
my $nb_of_groups = $arguments_tab[11];
my $include_exp_vec = $arguments_tab[12];
my $clustering_method = $arguments_tab[13];
my $sort_rows = $arguments_tab[14];
my $auto_scale = $arguments_tab[15];
my $zmin = $arguments_tab[16];
my $zmax = $arguments_tab[17];


print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 18)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/kmean_with_seqplots.R --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --bedVec $bed_vec --bedNameVec $bed_name_vec --organism $organism --genomeVersion $genome_version --binSize $bin_size --profileLengthBefore $profile_length_before --profileLengthAfter $profile_length_after --typeValue $type_value --outputFolder $output_folder --nbOfGroups $nb_of_groups --includeExpVec $include_exp_vec --clusteringMethod $clustering_method --sortRows $sort_rows --autoScale $auto_scale --zmin $zmin --zmax $zmax\n\n";

	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/kmean_with_seqplots.R --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --bedVec $bed_vec --bedNameVec $bed_name_vec --organism $organism --genomeVersion $genome_version --binSize $bin_size --profileLengthBefore $profile_length_before --profileLengthAfter $profile_length_after --typeValue $type_value --outputFolder $output_folder --nbOfGroups $nb_of_groups --includeExpVec $include_exp_vec --clusteringMethod $clustering_method --sortRows $sort_rows --autoScale $auto_scale --zmin $zmin --zmax $zmax";
	system($commandToLoad);
}else{
	
	# Define vector of colors
	my $col_value = $arguments_tab[18];
	
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/kmean_with_seqplots.R --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --bedVec $bed_vec --bedNameVec $bed_name_vec --organism $organism --genomeVersion $genome_version --binSize $bin_size --profileLengthBefore $profile_length_before --profileLengthAfter $profile_length_after --typeValue $type_value --outputFolder $output_folder --nbOfGroups $nb_of_groups --includeExpVec $include_exp_vec --clusteringMethod $clustering_method --sortRows $sort_rows --autoScale $auto_scale --zmin $zmin --zmax $zmax --colValue $col_value\n\n";

	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/kmean_with_seqplots.R --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --bedVec $bed_vec --bedNameVec $bed_name_vec --organism $organism --genomeVersion $genome_version --binSize $bin_size --profileLengthBefore $profile_length_before --profileLengthAfter $profile_length_after --typeValue $type_value --outputFolder $output_folder --nbOfGroups $nb_of_groups --includeExpVec $include_exp_vec --clusteringMethod $clustering_method --sortRows $sort_rows --autoScale $auto_scale --zmin $zmin --zmax $zmax --colValue $col_value";
	system($commandToLoad);
	
}





