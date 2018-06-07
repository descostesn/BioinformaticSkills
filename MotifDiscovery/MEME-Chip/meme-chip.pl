#!/usr/bin/env perl
#This script looks for motifs with meme-chip
# It uses the input file X.conf of format: outputfolder;motif_db;markovorder;thresholdclustering;thresholdtoincludemotif;analysisname;minimumotifswidth;maxmotifwidth;maxnbmotiftofind;minnbsitespermotif;maxnbsitepermotif;maxdatasize;dremeeval;maxnbmotiftofinddreme;minallowedmatchscore;maxregionsizecentrimo;thresholdevalcentrimo;fastafile 
#Nicolas Descostes dec 2016

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 18)
{
	die("Missing arguments for meme-chip.pl\n\n File should contain: outputfolder;motif_db;markovorder;thresholdclustering;thresholdtoincludemotif;analysisname;
	minimumotifswidth;maxmotifwidth;maxnbmotiftofind;minnbsitespermotif;maxnbsitepermotif;maxdatasize;dremeeval;maxnbmotiftofinddreme;minallowedmatchscore;maxregionsizecentrimo;thresholdevalcentrimo;fastafile\n\n");
}



my $output_folder= $arguments_tab[0]; 
my $motif_db= $arguments_tab[1];
my $markov_order= $arguments_tab[2];
my $threshold_clustering= $arguments_tab[3];
my $threshold_to_include_motif= $arguments_tab[4];
my $analysis_name= $arguments_tab[5];
my $minimum_motif_width= $arguments_tab[6];
my $maximum_motif_width= $arguments_tab[7];
my $maximum_number_motifs_to_find= $arguments_tab[8];
my $minimum_nb_sites_per_motif= $arguments_tab[9];     #This is the total number of sites in the training set where a single motif occurs
my $maximum_number_sites_per_motif= $arguments_tab[10];
my $max_dataset_size= $arguments_tab[11];
my $dreme_eval_threshold= $arguments_tab[12];
my $maximum_number_motifs_to_find_dreme= $arguments_tab[13];
my $minimum_allowed_match_score_centrimo= $arguments_tab[14];
my $max_region_size_centrimo= $arguments_tab[15];
my $threshold_eval_centrimo= $arguments_tab[16];
my $fasta_file= $arguments_tab[17];



print "This is job number $file_lineNumber \n";

print "meme-chip -o $output_folder -db $motif_db -dna -order $markov_order -norand -ccut 0 -group-thresh $threshold_clustering -filter-thresh $threshold_to_include_motif -desc $analysis_name -meme-minw $minimum_motif_width -meme-maxw $maximum_motif_width -meme-nmotifs $maximum_number_motifs_to_find -meme-minsites $minimum_nb_sites_per_motif -meme-maxsites $maximum_number_sites_per_motif -meme-p $nbCPU -meme-maxsize $max_dataset_size -dreme-e $dreme_eval_threshold -dreme-m $maximum_number_motifs_to_find_dreme -centrimo-local -centrimo-score  $minimum_allowed_match_score_centrimo -centrimo-maxreg $max_region_size_centrimo -centrimo-ethresh $threshold_eval_centrimo -centrimo-flip $fasta_file\n\n";


my $commandToLoad = "meme-chip -o $output_folder -db $motif_db -dna -order $markov_order -norand -ccut 0 -group-thresh $threshold_clustering -filter-thresh $threshold_to_include_motif -desc $analysis_name -meme-minw $minimum_motif_width -meme-maxw $maximum_motif_width -meme-nmotifs $maximum_number_motifs_to_find -meme-minsites $minimum_nb_sites_per_motif -meme-maxsites $maximum_number_sites_per_motif -meme-p $nbCPU -meme-maxsize $max_dataset_size -dreme-e $dreme_eval_threshold -dreme-m $maximum_number_motifs_to_find_dreme -centrimo-local -centrimo-score  $minimum_allowed_match_score_centrimo -centrimo-maxreg $max_region_size_centrimo -centrimo-ethresh $threshold_eval_centrimo -centrimo-flip $fasta_file";
system($commandToLoad);

