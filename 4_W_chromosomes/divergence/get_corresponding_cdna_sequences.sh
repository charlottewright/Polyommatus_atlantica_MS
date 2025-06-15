#!/usr/bin/env bash

# Read the input file
input_file="$1"

# Extract the first parts of the two lines starting with '>'
ID1=$(grep -m 1 '^>' "$input_file" | awk '{print $1}' | sed 's/^>//')
ID2=$(grep -m 2 '^>' "$input_file" | tail -n 1 | awk '{print $1}' | sed 's/^>//')

# Check if files matching the wildcard pattern exist
if ! ls cdna/"${ID1}"* >/dev/null 2>&1 || ! ls cdna/"${ID2}"* >/dev/null 2>&1; then
    echo "Error: One or both of the ID files do not exist"
    exit 1
fi

# Concatenate the matching files into a new file
output_file="${ID1}.plus_hit.fna"

cat nucleotide_sequences/"${ID1}"* nucleotide_sequences/"${ID2}"* > "$output_file"

# Provide feedback to the user
if [ $? -eq 0 ]; then
    echo "Files concatenated successfully into $output_file"
else
    echo "Error: Could not concatenate files"
    exit 1
fi
