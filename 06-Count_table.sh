#!/bin/bash

create an array with the names of all sample files
samples=(*.counts.txt)

# Initialize a variable to hold the output file name
output_file="Plutella-mapped_counts.txt"

# Loop through each sample and extract the 1st and 7th columns using awk
for sample in "${samples[@]}"
do
    if [ -f "${output_file}" ]; then
        # Append data to the output file if it already exists
        paste ${output_file} <(awk 'BEGIN {OFS="\t"} {print $7}' "${sample}") | grep -v '^\#' > temp.txt
        mv temp.txt ${output_file}
    else
        # Create the output file with data from the first sample
        awk 'BEGIN {OFS="\t"} {print $1,$7}' "${sample}" > ${output_file}
    fi
done

echo "Done!"









