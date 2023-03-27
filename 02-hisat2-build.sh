
#!/bin/bash

#SBATCH --job-name=hisat2 # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks-per-node=16 # Number of CPU cores per nod
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written
#SBATCH --mail-user=nasir.abbasi@univ-tours.fr
#SBATCH --mail-type=begin
#SBATCH --mail-type=end




#module load hisat2

GENOME_FNA="/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Plutella_reference/genome_assemblies_genome_fasta/ncbi-genomes-2022-08-03/GCA_019096205.1_PxLV.1_genomic.fna/GCA_019096205.1_PxLV.1_genomic.fna"
GENOME=${GENOME_FNA%.*}                # Drops the .fna extension
hisat2-build ${GENOME_FNA} ${GENOME}
