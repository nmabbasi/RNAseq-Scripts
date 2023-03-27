#!/bin/bash

#SBATCH --job-name=trimmomatic # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --cpus-per-task=16 # number of cpus to use
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --array=1-1
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written



#activate required conda env for fastp
conda activate TRIMMOMATIC

# set up input and output directories
input_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc"
output_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc/trim_output"


#Get the input file names for the current array task
input_files=($(ls ${input_dir}/*_1.fq.gz))
file_name=${input_files[$SLURM_ARRAY_TASK_ID-1]}
base_name=$(basename ${file_name} _1.fq.gz)

        # Run trimmomatic to trim the reads
trimmomatic PE \
  -threads ${SLURM_CPUS_PER_TASK} \
  -phred33 \
  ${input_dir}/${base_name}_1.fq.gz \
  ${input_dir}/${base_name}_2.fq.gz \
  ${output_dir}/${base_name}_1_trimmed.fq.gz \
  ${output_dir}/${base_name}_1_unpaired.fq.gz \
  ${output_dir}/${base_name}_2_trimmed.fq.gz \
  ${output_dir}/${base_name}_2_unpaired.fq.gz \
  ILLUMINACLIP:/home/biopatic/SOFTS/Miniconda/envs/TRIMMOMATIC/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36 \
  AVGQUAL:25 \
  HEADCROP:10 \
  CROP:90 \
  AVGQUAL:20 \
  MAXINFO:90:0.1 \
  TOPHRED33


echo "Finished trimming ${base_name} at $(date)"

