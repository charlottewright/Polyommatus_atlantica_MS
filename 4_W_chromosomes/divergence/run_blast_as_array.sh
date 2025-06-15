#!/bin/bash

# running in an array, using $LSB_JOBINDEX
name=$(sed -n ${LSB_JOBINDEX}p list_of_input_files.txt)

# run blastp
echo ${name}
blastp -query ${name} -db ../proteins_except_Ws -evalue 0.00001 -outfmt 6 -out ${name}_W_blast_matches.txt
