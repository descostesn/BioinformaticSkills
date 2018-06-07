###############
# This script uses the bioconductor package DiffBind to perform differential 
# binding analysis with ChIP-Seq data. This script is written specifically and
# limited to the study of Gary, i.e, it handles only two conditions and the 
# contrasts are built on these two.
# conf file should contain: CSVFILE;ANALYSISNAME;ELONGATIONSIZE;OUTPUTFOLDER;SCALINGFACTORS
# Descostes MAY 2018
###############

use strict;
my $file_lineNumber = $ARGV[0];
my $input_file = $ARGV[1];


#Retrieve the parameters
my $file_path = `head -n $file_lineNumber $input_file | tail -n1`;
chomp $file_path;

my @arguments_tab = split(';', $file_path);

if(scalar(@arguments_tab) != 5)
{
    die("Missing arguments for diffbind.pl\n\n File should contain: CSVFILE;ANALYSISNAME;ELOLNGATIONSIZE;OUTPUTFOLDER;SCALINGFACTORS\n\n");
}

my $csv_file = $arguments_tab[0];
my $analysis_name = $arguments_tab[1];
my $elongation_size = $arguments_tab[2];
my $output_folder = $arguments_tab[3];
my $scaling_factors = $arguments_tab[4];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/diffbind.R --csvFile $csv_file --analysisName $analysis_name --elongationSize $elongation_size --outputFolder $output_folder --scalingFactors $scaling_factors\n\n";
my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/diffbind.R --csvFile $csv_file --analysisName $analysis_name --elongationSize $elongation_size --outputFolder $output_folder --scalingFactors $scaling_factors";
system($commandToLoad);
