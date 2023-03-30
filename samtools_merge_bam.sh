#!/bin/bash

# options for sbatch
#SBATCH --job-name=merge_bam_index # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # number of cpus to use
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
#SBATCH --array=1-24
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written

set -o xtrace

# Set input and output directories
in_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Enterobacter_reference/Enterobacter/3-Enterobacter_Mapping/sorted_bam"
out_dir="/home/biopatic/abbasi/P4--RNA-Seq-loreta/Ref_Database/Loretta/Enterobacter_reference/Enterobacter/3-Enterobacter_Mapping/sorted_bam/merge_sorted_and_index"

#Define sample names and corresponding input files
SAMPLE_NAMES=("I2" "I4" "I6" "I8" "M3" "M4" "M6" "M8" "R3" "R7" "R8" "R9" "S2" "S3" "S4" "S5" "T3" "T4" "T6" "T8" "X1" "X3" "X4" "X7")



for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
	OUTPUT_FILE="${out_dir}/${SAMPLE_NAME}_sorted.bam"
    

	if [ -f "${OUTPUT_FILE}" ]; then
        	echo "${OUTPUT_FILE} already exists. Skipping merge and index for ${SAMPLE_NAME}."
	else
	    INPUT_FILES=()
     	    for LANE_NUMBER in {1..3}; do
        		# Search for files matching the pattern
    
			FILE_NAME="${SAMPLE_NAME}_FKRN*_*_L${LANE_NUMBER}_sorted.bam"
        	while read -r -d '' FULL_FILE_NAME; do
            		INPUT_FILES+=("${FULL_FILE_NAME}")
        	done < <(find "${in_dir}" -name "${FILE_NAME}" -type f -print0)
    	    done






    		# Merge sorted BAM files using samtools
    		OUTPUT_FILE="${out_dir}/${SAMPLE_NAME}_sorted.bam"
    		samtools merge -@ 16 "${OUTPUT_FILE}" "${INPUT_FILES[@]}"

    		# Index the merged sorted BAM file
    		samtools index "${OUTPUT_FILE}"
	fi
done
