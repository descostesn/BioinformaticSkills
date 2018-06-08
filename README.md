# EMBL-Showcase

## Introduction


This repository aims at explaining how my work is organized and to provide a list of tools that I use for analyzing sequencing data. By clicking names, one can visualize the corresponding code.

The folder [ChIPSeqSpike_dev](https://github.com/descostesn/embl-showcase/tree/master/ChIPSeqSpike_dev) provides the current code under development for [ChIPSeqSpike](https://www.bioconductor.org/packages/devel/bioc/html/ChIPSeqSpike.html). This code is the most up-to-date in terms of my R knowledge. Other scripts provided are often developed for internal use and do not intend to use the optimal solutions necessarily. The provided code was not 'polished' and is shown as being developed on a daily basis. As an embeded bioinformatician, I always try to find the best compromise between code optimization and time of development.

## General Overview

### System spec

All codes were developed on a local computer under **Ubuntu operating system**. I use **Eclipse** Oxygen 3a Release 4.7.3a committing with **Egit** to **Github**. I make use of specific perspectives (StatET, Perl, Java, Bash) in Eclipse to **optimize development** and I am using code templates as well to increase my productivity.</p>
The **HPC** system used (via scp, ssh) is [Phoenix](https://genome.med.nyu.edu/hpcf/wiki/Manual:Cluster_User_Guide#Foreword), CentOS release 6, that makes use of **SGE**. It is composed of 2 head nodes and 70 compute nodes. 7 are equipped with a total of 21 GPUs, and one is a high-memory node (1 TB of RAM).

### Code and Jobs

I am making use of R and Bioconductor extensively. Most of my scripts are structured in the following manner:

  1. A brief description and the date (MM/YYYY) of development
  2. A parameter section with a description of each argument used
  3. A functions section
  4. A main section

To process data in batch, I am usually building '.conf' files containing all parameters separated by semicolons (see [code generation repository](https://github.com/descostesn/embl-showcase/tree/master/Code_generation) for few examples). To feed script with parameters, I am using perl scripts that call R script through the `system(Rscript path*.R)` command. The perl script itself is called by a bash script defining hpc parameters. This latter is qsub to the grid. In summary, my data treatments are organized as follows:

  1. Bash defining job scheduler parameters taking a *.conf* file as input
  2. Perl script formating parameters and calling the *processing* script
  3. Processing script


## Sequencing data processing and analysis

**Click on the process that you are interested in to see the code:**

  1. **For the latest and most representative code:** [ChipSeqSpike development](https://github.com/descostesn/embl-showcase/tree/master/ChIPSeqSpike_dev)
  2. A [pipeline](https://github.com/descostesn/embl-showcase/tree/master/fastq2Bigwigs) for ChIP-seq and RNA-Seq 
  3. Demultiplexing and trimming: with [BB](https://github.com/descostesn/embl-showcase/tree/master/Demultiplexing_trimming/BBprimers), [Illumina primers](https://github.com/descostesn/embl-showcase/tree/master/Demultiplexing_trimming/IlluminaPrimers), [VLB](https://github.com/descostesn/embl-showcase/tree/master/Demultiplexing_trimming/VLB), [cutadapt](https://github.com/descostesn/embl-showcase/tree/master/Demultiplexing_trimming/cutadapt)
  4. QC control and filtering: [fastqc](https://github.com/descostesn/embl-showcase/tree/master/QualityControlFiltering), [NGSQCToolkit](https://github.com/descostesn/embl-showcase/tree/master/QualityControlFiltering), 
  5. Alignment: [Bowtie](https://github.com/descostesn/embl-showcase/tree/master/Alignment/Bowtie), [Bowtie2](https://github.com/descostesn/embl-showcase/tree/master/Alignment/Bowtie2), [STAR](https://github.com/descostesn/embl-showcase/tree/master/Alignment/STAR), [Tophat](https://github.com/descostesn/embl-showcase/tree/master/Alignment/Tophat)
  6. Bam operations: [ALL](https://github.com/descostesn/embl-showcase/tree/master/BamOperations), [bamutils](https://github.com/descostesn/embl-showcase/tree/master/BamOperations/bamutils), [samtools](https://github.com/descostesn/embl-showcase/tree/master/BamOperations/samtools), [picardtools](https://github.com/descostesn/embl-showcase/tree/master/BamOperations/picardtools) 
  7. Pasha Pipeline (published in Bioinformatics journal: [here](https://github.com/descostesn/embl-showcase/tree/master/Pasha))
  8. Format conversion (not exhaustive): [sratoolkit](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/sratoolkit), [bedtools](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/bedtools), [samtools](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/samtools), [bedops](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/bedops), [bed2GFF](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/), [wigVar2Fix](https://github.com/descostesn/embl-showcase/tree/master/FormatConversion/)
  9. Visualization: [Boxplots](https://github.com/descostesn/embl-showcase/tree/master/Visualization/Boxplots),[UCSC tracks](https://github.com/descostesn/embl-showcase/tree/master/Visualization/UCSCtracks), [linear regression](https://github.com/descostesn/embl-showcase/tree/master/Visualization/LinearRegression), [MAPlot](https://github.com/descostesn/embl-showcase/tree/master/Visualization/MAplot), [PCA](https://github.com/descostesn/embl-showcase/tree/master/Visualization/PCA), [profiles](https://github.com/descostesn/embl-showcase/tree/master/Visualization/Profiling), [venn diagrams](https://github.com/descostesn/embl-showcase/tree/master/Visualization/VennDiagrams)
  10. Interval operations: [Centering](https://github.com/descostesn/embl-showcase/tree/master/IntervalOperations/Centering), [union and intersection](https://github.com/descostesn/embl-showcase/tree/master/IntervalOperations/IntersectionUnion)
  11. Clustering: [seqplots](https://github.com/descostesn/embl-showcase/tree/master/Clustering/seqplots), [deeptools](https://github.com/descostesn/embl-showcase/tree/master/Clustering/deeptools)
  12. Differential binding: [diffbind](https://github.com/descostesn/embl-showcase/tree/master/Differential/ChIP-Seq/DiffBind), [macs2diff](https://github.com/descostesn/embl-showcase/tree/master/Differential/ChIP-Seq/Macs2Diff), [MAnorm](https://github.com/descostesn/embl-showcase/tree/master/Differential/ChIP-Seq/MAnorm)
  13. Gene Ontologies: [chipEnrich](https://github.com/descostesn/embl-showcase/tree/master/GeneOntologies/ChIPEnrich), [clusterProfiler](https://github.com/descostesn/embl-showcase/tree/master/GeneOntologies/clusterProfiler), [Gage](https://github.com/descostesn/embl-showcase/tree/master/GeneOntologies/Gage)
  14. Markov model: [chromHMM](https://github.com/descostesn/embl-showcase/tree/master/MarkovModel/chromHMM), [Model Analysis](https://github.com/descostesn/embl-showcase/tree/master/MarkovModel/chromHMM_analysis)
  15. Motif discovery: [g-quad](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/Gquads), [gimmemotifs](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/gimmemotifs), [meme](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/MEME), [meme-chip](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/MEME-Chip), [seqpattern](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/seqpattern-modified), [rsat](https://github.com/descostesn/embl-showcase/tree/master/MotifDiscovery/rsat)
  16. Peak calling: [GEM](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/GEM), [hiddenDomains](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/hiddenDomains), [macs2](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/macs2), [PeakSeq](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/PeakSeq), [SICER](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/SICER), [SPP](https://github.com/descostesn/embl-showcase/tree/master/PeakCalling/SPP)
  17. RNA-Seq: [DESeq2](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/DESeq2), [edgeR](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/edgeR), [maSigPro](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/maSigPro), [TCSeq](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/TCSeq), [DEXSeq](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/DEXSeq), [Tissue specificity](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/TissueSpecificity), [interval division](https://github.com/descostesn/embl-showcase/tree/master/Differential/RNASeq/DivisionCategory)

## Java interface for Biologists

[Genomic Box](https://github.com/descostesn/embl-showcase/tree/master/JavaGUIForBiologists)
