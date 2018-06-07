# EMBL-showcase

## Introduction


This repository aims at explaining how my work is organized and to provide a list of tools that I use for analyzing sequencing data.

<p style="text-align: justify;"> The folder 'bioconductor-development' provides the current code under development for [ChIPSeqSpike](https://www.bioconductor.org/packages/devel/bioc/html/ChIPSeqSpike.html). This code is the most up-to-date in terms of my R knowledge. Other scripts provided are often developed for internal use and do not intend to use the optimal solutions necessarily. As an embeded bioinformatician, I always try to find the best compromise between code optimization and time of development.</p>

## General Overview

### System spec
<p style='text-align: justify;'> All the code has been developed on a local computer under Ubuntu operating system. I use Eclipse Oxygen 3a Release 4.7.3a committing with Egit to Github. I make use of specific perspectives (StatET, Perl, Java, Bash) in Eclipse to optimize development and I am using code templates as well to increase my productivity.</p>
The HPC system used (via scp, ssh) is [Phoenix](https://genome.med.nyu.edu/hpcf/wiki/Manual:Cluster_User_Guide#Foreword), CentOS release 6, that makes use of SGE. It is composed of 2 head nodes and 70 compute nodes. 7 are equipped with a total of 21 GPUs, and one is a high-memory node (1 TB of RAM).</p>

### Code and Jobs

I am making use of R and Bioconductor extensively. All my scripts are structured in the following manner:

  1. A brief description and the date (MM/YYYY) of development. 
  2. A parameter section with a description of each argument used
  3. A functions section
  4. A main section

<p style='text-align: justify;'> To process data in batch, I am usually building '.conf' files containing all parameters separated by semicolons (see [code generation repository]()). To feed script with parameters, I am using perl scripts that call R script through the `system(Rscript path*.R)` command.</p>
