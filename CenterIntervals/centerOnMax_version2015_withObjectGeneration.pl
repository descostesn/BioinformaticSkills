###############
# This script combines the object generation (with save) and centering on max
# conf file should contain: see script
# Descostes March 2018
###############

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
    die("Missing arguments for centerOnMax_version2015_withObjectGeneration.pl\n\n File should contain: see script\n\n");
}


my $annotationFolder = $arguments_tab[0];
my $annotationFile = $arguments_tab[1];
my $annotationComment = $arguments_tab[2];
my $annotationFormat = $arguments_tab[3];
my $wigFolders = $arguments_tab[4];
my $wigFiles = $arguments_tab[5];
my $expNames = $arguments_tab[6];
my $annotationExtensions = $arguments_tab[7];
my $spaceSizes = $arguments_tab[8];
my $outputFolders = $arguments_tab[9];
my $binSize = $arguments_tab[10];


print "This is job number $file_lineNumber \n";

print "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/centerOnMax_version2015_withObjectGeneration.R --annotationFolder $annotationFolder --annotationFile $annotationFile --annotationComment $annotationComment --annotationFormat $annotationFormat --wigFolders $wigFolders --wigFiles $wigFiles --expNames $expNames --annotationExtensions $annotationExtensions --spaceSizes $spaceSizes --outputFolders $outputFolders --binSize $binSize --nbCPU $nbCPU";

my $commandToLoad = "Rscript /ifs/home/descon01/cluster/scripts/R_scripts/centerOnMax_version2015_withObjectGeneration.R --annotationFolder $annotationFolder --annotationFile $annotationFile --annotationComment $annotationComment --annotationFormat $annotationFormat --wigFolders $wigFolders --wigFiles $wigFiles --expNames $expNames --annotationExtensions $annotationExtensions --spaceSizes $spaceSizes --outputFolders $outputFolders --binSize $binSize --nbCPU $nbCPU";
system($commandToLoad);

