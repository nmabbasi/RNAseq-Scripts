/*
RNAseq pipeline
*/

// Specify the version of Nextflow
nextflow.version = '23.04.1.5866'

params {
    path_fastq = "data/*_{1,2}.fastq.gz"
    ref_genome = "genome.fa"
    gtf_file = "annotations.gtf"
}

// Define the process for trimming
process trim_reads {
    input:
    file fastq from path_fastq

    output:
    file "trimmed_${fastq.baseName}" into trimmed_fastqs

    script:
    """
    fastp --in1 ${fastq} --in2 ${fastq.baseName.replace("1.fastq.gz", "2.fastq.gz")} \
    --out1 trimmed_${fastq.baseName} --out2 trimmed_${fastq.baseName.baseName.replace("1.fastq.gz", "2.fastq.gz")} \
    --thread 4 --html trimmed_${fastq.baseName}.html --json trimmed_${fastq.baseName}.json
    """
}

// Define the process for mapping with HISAT2
process hisat2_mapping {
    input:
    file fastq from trimmed_fastqs
    path genome from ref_genome

    output:
    file "aligned_${fastq.baseName}.bam" into mapped_bams

    script:
    """
    hisat2 -p 4 -x ${genome} -1 ${fastq.baseName} -2 ${fastq.baseName.replace("1.fastq.gz", "2.fastq.gz")} \
    | samtools view -b - > aligned_${fastq.baseName}.bam
    """
}

// Define the process for sorting and indexing with Samtools
process sort_and_index {
    input:
    file bam from mapped_bams

    output:
    file "sorted_${bam.baseName}.bam" into sorted_bams

    script:
    """
    samtools sort -o sorted_${bam.baseName}.bam ${bam}
    samtools index sorted_${bam.baseName}.bam
    """
}

// Define the process for quantification with FeatureCounts
process feature_counts {
    input:
    file bam from sorted_bams
    path annotations from gtf_file

    output:
    file "counts_${bam.baseName}.txt" into count_files

    script:
    """
    featureCounts -T 4 -a ${annotations} -o counts_${bam.baseName}.txt ${bam}
    """
}

// Define the process for performing Deseq2 analysis with R
process deseq2_analysis {
    input:
    file counts from count_files

    output:
    file "deseq2_results.txt" into results

    script:
    """
    library(DESeq2)
    counts <- read.table("${counts}", header=TRUE, row.names=1, sep="\t")
    coldata <- data.frame(condition=c(rep("treated",3), rep("untreated",3)))
    rownames(coldata) <- colnames(counts)
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    write.table(resOrdered, file="deseq2_results.txt", sep="\t")
    """
}

// Define the pipeline
workflow {
    // Define the pipeline inputs
    input_fastq = set path_fastq

    // Trim reads
    trimmed_fastqs = trim_reads(input_fastq)

    // Map reads with HISAT2
    mapped_bams

    // sorting and indexing with Samtools
    sort_and_index

   // quantification with FeatureCounts
   feature_counts

  // performing Deseq2 analysis with R
  deseq2_analysis

