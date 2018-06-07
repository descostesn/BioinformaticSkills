#!/usr/bin/env perl
#This script generate heatmaps and av profiles of motifs from a fasta file.
# It uses the input file files_wig_X.conf of format: fastafile;flankup;flankdown;nbin;outputfolder;consensusiupacvec;consensustoolsvec;consensusmotifvec;pwmvec;pwmnamevec;minscorepwm
#Nicolas Descostes january 2017

use strict;
my $file_lineNumber = $ARGV[0];
my $nbCPU = $ARGV[1];
my $input_file = $ARGV[2];

#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 11)
{
	die("Missing arguments for seqpattern.pl\n\n File should contain: fastafile;flankup;flankdown;nbin;outputfolder;consensusiupacvec;consensustoolsvec;consensusmotifvec;pwmvec;pwmnamevec;minscorepwm\n\n");
}


my $fasta_file = $arguments_tab[0];
my $flank_up = $arguments_tab[1];
my $flank_down = $arguments_tab[2];
my $n_bin = $arguments_tab[3];
my $output_folder = $arguments_tab[4];
my $consensus_iupac_vec = $arguments_tab[5];
my $consensus_tools_vec = $arguments_tab[6];
my $consensus_motif_name = $arguments_tab[7];
my $pwm_vec = $arguments_tab[8];
my $pwm_name_vec = $arguments_tab[9];
my $min_score_pwm = $arguments_tab[10];



print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/seqpattern.R --fastaFile $fasta_file --flankUp $flank_up --flankDown $flank_down --Nbin $n_bin --outputFolder $output_folder --nbCpu $nbCPU --consensusIupacVec $consensus_iupac_vec --consensusToolsVec $consensus_tools_vec --consensusMotifName_vec $consensus_motif_name --pwmVec $pwm_vec --pwmNameVec $pwm_name_vec --minScorePwm $min_score_pwm\n\n";


my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/seqpattern.R --fastaFile $fasta_file --flankUp $flank_up --flankDown $flank_down --Nbin $n_bin --outputFolder $output_folder --nbCpu $nbCPU --consensusIupacVec $consensus_iupac_vec --consensusToolsVec $consensus_tools_vec --consensusMotifName_vec $consensus_motif_name --pwmVec $pwm_vec --pwmNameVec $pwm_name_vec --minScorePwm $min_score_pwm";
system($commandToLoad);



 


