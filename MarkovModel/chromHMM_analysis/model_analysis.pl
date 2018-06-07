###############
# This script performs different analysis and statistics on the defined markov model.
# conf file should contain: BEDFILEMODEL;OUTPUTFOLDER;GENOMEVERSION;REFSEQFILE;BAMFILERNASEQ;SINGLEEND;OUTPUTFORMAT;COLVECORIGINAL
# Descostes June 2017
###############

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 9)
{
	die("Missing arguments for model_analysis.pl\n\n File should contain: BEDFILEMODEL;OUTPUTFOLDER;GENOMEVERSION;REFSEQFILE;BAMFILERNASEQ;SINGLEEND;OUTPUTFORMAT;COLVECORIGINAL\n\n");
}

my $bed_file_model = $arguments_tab[0];
my $output_folder = $arguments_tab[1];
my $genome_version = $arguments_tab[2];
my $refseq_file = $arguments_tab[3];
my $refseq_anno_gtf = $arguments_tab[4];
my $bam_file_rnaseq = $arguments_tab[5];
my $single_end = $arguments_tab[6];
my $output_format = $arguments_tab[7];
my $col_vec_original = $arguments_tab[8];

print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/model_analysis.R --bedFileModel $bed_file_model --outputFolder $output_folder --genomeVersion $genome_version --refseqFile $refseq_file --refseqAnnoGtf $refseq_anno_gtf --bamFileRnaseq $bam_file_rnaseq --singleEnd $single_end --outputFormat $output_format --colVecOriginal $col_vec_original\n\n";
my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/model_analysis.R --bedFileModel $bed_file_model --outputFolder $output_folder --genomeVersion $genome_version --refseqFile $refseq_file --refseqAnnoGtf $refseq_anno_gtf --bamFileRnaseq $bam_file_rnaseq --singleEnd $single_end --outputFormat $output_format --colVecOriginal $col_vec_original";
system($commandToLoad);






