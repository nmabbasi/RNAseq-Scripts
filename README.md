# RNAseq-Scripts
This is RNAseq scripts which I used for the analysis of 24 samples. 

The samples were processed in different lanes. so, I created a script and I merged bam files later on before applying FeatureCounts.

I did following steps:

Extracting data online using wget on HPC slurm cluster
01. Quality check: fastqc
02. Indexing reference:hisat2-build
03. Trimming using fastp and Trimmomatic
04. Mapping using Hisat2
05. FeatureCounts
06. Creating a count-table

