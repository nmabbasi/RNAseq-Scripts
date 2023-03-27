#!/bin/bash

#SBATCH --job-name=hisat2 # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16 # Number of CPU cores per nod
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --array=2-46
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written
#SBATCH --mail-user=nasir.abbasi@univ-tours.fr
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -o xtrace


# Path to reference genome and index prefix
REF_GENOME=/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/2-Plutella_genome/GCA_019096205.1_PxLV.1_genomic.fna
INDEX=/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/2-Plutella_index/GCA_019096205.1_PxLV.1_genomic

# Path to input fastq files and output directories
input_dir=/home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc/2-fastp_Trimming_output/fastqc_of_trimmedSamples
OUTDIR=/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/3-Plutella_Mapping/BAM_files

sorted_bam=/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/3-Plutella_Mapping/sorted_bam

#Get the input file names for the current array task
input_files=($(ls ${input_dir}/*_1_trimmed.fq.gz))
file_name=${input_files[$SLURM_ARRAY_TASK_ID-1]}
base_name=$(basename ${file_name} _1_trimmed.fq.gz)

# Run HISAT2
hisat2 -x $INDEX -1 $input_dir/${base_name}_1_trimmed.fq.gz -2 ${input_dir}/${base_name}_2_trimmed.fq.gz -S $OUTDIR/${base_name}.sam -p 16 --un-conc-gz $OUTDIR/unmapped_${base_name} 
   

	 # Convert SAM to BAM 
    	samtools view -bS $OUTDIR/${base_name}.sam >$OUTDIR/${base_name}.bam

	#Sort BAM file
	 samtools sort  $OUTDIR/${base_name}.bam $sorted_bam/${base_name}_sorted 

    	# Index BAM file
    	samtools index $sorted_bam/${base_name}_sorted.bam 

    

