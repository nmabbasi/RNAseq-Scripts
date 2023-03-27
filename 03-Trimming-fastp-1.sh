#!/bin/bash

# options for sbatch
#SBATCH --job-name=fastp # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --cpus-per-task=16 # number of cpus to use
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --array=1-1
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written



#activate required conda env for fastp
conda activate fastp

# set up input and output directories
input_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc"
output_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/X204SC22123284-Z01-F001/01.RawData/1-All_samples_fastqc/fastp_Test"


#Get the input file names for the current array task
input_files=($(ls ${input_dir}/*_1.fq.gz))
file_name=${input_files[$SLURM_ARRAY_TASK_ID-1]}
base_name=$(basename ${file_name} _1.fq.gz)

    	# Run Fastp to trim the reads
fastp -i ${input_dir}/${base_name}_1.fq.gz \
      -I ${input_dir}/${base_name}_2.fq.gz \
      -o ${output_dir}/${base_name}_1_trimmed1.fq.gz \
      -O ${output_dir}/${base_name}_2_trimmed1.fq.gz \
	--detect_adapter_for_pe --qualified_quality_phred 20 --average_qual 20 --low_complexity_filter  --overrepresentation_analysis --dedup --correction --thread 16 --json ${output_dir}/${base_name}.json --html ${output_dir}/${base_name}.html 	

echo "Finished trimming ${base_name} at $(date)"

