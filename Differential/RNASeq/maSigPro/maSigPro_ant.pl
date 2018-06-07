###############
# This script aims at performing time course rnaseq analysis.
# conf file should contain: BAMFILEVEC;REFSEQANNO;EXPNAMESVEC;DESIGNFILE;OUTPUTFOLDER;MINCLUSTERNB
# Descostes May 2017
###############

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 6)
{
	die("Missing arguments for maSigPro_ant.pl\n\n File should contain: BAMFILEVEC;REFSEQANNO;EXPNAMESVEC;DESIGNFILE;OUTPUTFOLDER;MINCLUSTERNB\n\n");
}

my $bam_files_vec   = $arguments_tab[0];
my $refseq_anno  = $arguments_tab[1];
my $expnames_vec  = $arguments_tab[2];
my $design_file  = $arguments_tab[3];
my $output_folder  = $arguments_tab[4];
my $min_cluster_nb  = $arguments_tab[5];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maSigPro_ant.R --bamFilesVec $bam_files_vec --refseqAnno $refseq_anno --expnamesVec $expnames_vec --designFile $design_file --outputFolder $output_folder --minClusterNb $min_cluster_nb\n\n";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/maSigPro_ant.R --bamFilesVec $bam_files_vec --refseqAnno $refseq_anno --expnamesVec $expnames_vec --designFile $design_file --outputFolder $output_folder --minClusterNb $min_cluster_nb";
system($commandToLoad);









