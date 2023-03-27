#!/bin/bash

#SBATCH --job-name=featureCounts # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16 # Number of CPU cores per nod
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written
#SBATCH --mail-user=nasir.abbasi@univ-tours.fr
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


set -o xtrace

cd /home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA 

ANNOT_GFF="/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/3-Plutella_Annotation/GCA_019096205.1_PxLV.1_genomic.gff"          # file should be decompressed, e.g., gzip -d {filename}.gff.gz
BAM_DIR="/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/2-Plutella_Mapping/sorted_bam/merged_index_sorted_file"
OUT_DIR=/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/PLUTELLA/4-counts/count-by-gene-ID

for bam_file in ${BAM_DIR}/*_sorted.bam
do
# Get the sample name from the bam file name
    sample_name=$(basename ${bam_file})

	#Run Feature count
	
	featureCounts -T 16 -p -t gene -B -C -g ID \
	-F GFF -a ${ANNOT_GFF} \
	-o ${OUT_DIR}/${sample_name}.counts.txt \
	${bam_file}
	
done

