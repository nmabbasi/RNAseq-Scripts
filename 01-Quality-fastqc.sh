#!/bin/bash

#SBATCH --job-name=fastqc # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks-per-node=48 # Number of CPU cores per nod
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nasir.abbasi@univ-tours.fr




#load required modules
module load fastqc

# Change to directory with fastq files
#cd /home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc

# Run FastQC on all fastq files
for file in *.fq.gz; do
  fastqc -o /home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc/Output_all_samples "$file"
done
