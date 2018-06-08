###############
# This script aims at performing venn diagram of interval overlaps between gff files. The script takes from 2 to 4 files.
# conf file should contain: GFFFILEVEC;OUTPUTFOLDER;COMPARISONTITLE;EXPNAMETAB;GENOMEVERSION;UPSTREAMHEATMAP;DOWNSTREAMHEATMAP;BIGWIGVEC;BIGWIGNAMEVEC;ORGANISMNAME;CENTERCRITERION;COLVEC;OUTPUTFORMAT or if using max: GFFFILEVEC;OUTPUTFOLDER;COMPARISONTITLE;EXPNAMETAB;GENOMEVERSION;UPSTREAMHEATMAP;DOWNSTREAMHEATMAP;BIGWIGVEC;BIGWIGNAMEVEC;ORGANISMNAME;CENTERCRITERION;COLVEC;REFEXP
# Descostes feb 2016
###############

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 9 && scalar(@arguments_tab) != 13 && scalar(@arguments_tab) != 14)
{
	die("Missing arguments for vennDiagram_overlapGFF.pl\n\n File should contain: GFFFILEVEC;OUTPUTFOLDER;COMPARISONTITLE;EXPNAMETAB;GENOMEVERSION;ORGANISMNAME;CENTERCRITERION;COLVEC;OUTPUTFORMAT\n 
	or if using bigwig: GFFFILEVEC;OUTPUTFOLDER;COMPARISONTITLE;EXPNAMETAB;GENOMEVERSION;ORGANISMNAME;CENTERCRITERION;COLVEC;OUTPUTFORMAT;UPSTREAMHEATMAP;DOWNSTREAMHEATMAP;BIGWIGVEC;BIGWIGNAMEVEC\n
	or if using max: GFFFILEVEC;OUTPUTFOLDER;COMPARISONTITLE;EXPNAMETAB;GENOMEVERSION;ORGANISMNAME;CENTERCRITERION;COLVEC;OUTPUTFORMAT;UPSTREAMHEATMAP;DOWNSTREAMHEATMAP;BIGWIGVEC;BIGWIGNAMEVEC;REFEXP\n\n");
}

my $gff_file_vec = $arguments_tab[0];
my $output_folder = $arguments_tab[1];
my $comparison_title = $arguments_tab[2];
my $expnames_tab = $arguments_tab[3];
my $genome_version = $arguments_tab[4];
my $organism_name = $arguments_tab[5];
my $center_criterion = $arguments_tab[6];
my $col_vec = $arguments_tab[7];
my $output_format = $arguments_tab[8];

print "This is job number $file_lineNumber \n";

if(scalar(@arguments_tab) == 9)
{
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format";
	system($commandToLoad);

}elsif(scalar(@arguments_tab) == 13){

	my $upstream_heatmap = $arguments_tab[9];
	my $downstream_heatmap = $arguments_tab[10];
	my $bigwig_vec = $arguments_tab[11];
	my $bigwig_name_vec = $arguments_tab[12];

	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format --upstreamHeatmap $upstream_heatmap --downstreamHeatmap $downstream_heatmap --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format --upstreamHeatmap $upstream_heatmap --downstreamHeatmap $downstream_heatmap --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec";
	system($commandToLoad);
	 

}else{
	
	my $upstream_heatmap = $arguments_tab[9];
	my $downstream_heatmap = $arguments_tab[10];
	my $bigwig_vec = $arguments_tab[11];
	my $bigwig_name_vec = $arguments_tab[12];
	my $ref_exp_max_centering = $arguments_tab[13];
	
	print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format --upstreamHeatmap $upstream_heatmap --downstreamHeatmap $downstream_heatmap --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --refExpMaxCentering $ref_exp_max_centering\n\n";
	my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/vennDiagram_overlapGFF.R --gffFileVec $gff_file_vec --outputFolder $output_folder --comparisonTitle $comparison_title --expnamesTab $expnames_tab --genomeVersion $genome_version --organismName $organism_name --centerCriterion $center_criterion --colVec $col_vec --outputFormat $output_format --upstreamHeatmap $upstream_heatmap --downstreamHeatmap $downstream_heatmap --bigwigVec $bigwig_vec --bigwigNameVec $bigwig_name_vec --refExpMaxCentering $ref_exp_max_centering";
	system($commandToLoad);
}
