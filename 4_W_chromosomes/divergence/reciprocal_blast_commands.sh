#!/usr/bin/bash

# Aim: calculate divergence between W and Z portions

# convert from multi-line to single line per protein
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.fa > BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.fa

# filter for longest isoform per protein at the start of the process
mamba activate AGAT_env
agat_sp_keep_longest_isoform.pl -gff BRAKER_1_BRAKER_2_TSEBRA.output.gtf -o BRAKER_1_BRAKER_2_TSEBRA.output.filtered.gtf

# use filtered gtf to filter protein file
grep -A1 -f <(grep ID BRAKER_1_BRAKER_2_TSEBRA.output.filtered.gtf | grep 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq) BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.fa > BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa


# Strategy is to blast all Z proteins to everything else, then blast W to everything else
# ask which output hits have same TOP blast hit i.e. best hit i.e. intersection between them
# using this as output solves problem of getting multiple hits per protein.

# filter out Z-linked genes
sed -e '/SUPER_Z/,+1d' BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa > BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.noZs.fa
# filter out W-linked genes
sed -e '/SUPER_W/,+1d' BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa > BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.noWs.fa

# get just Z-linked genes
grep -A1 'SUPER_Z' BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa > Z_proteins.fa
# get just W-linked genes
grep -A1 'SUPER_W' BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa > W_proteins.fa

# make blast databases
mamba activate blast_env
makeblastdb -in BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.noZs.fa -dbtype prot -out proteins_except_Zs
makeblastdb -in BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.noWs.fa -dbtype prot -out proteins_except_Ws

# Step 2: Make each Z gene a query
mkdir Z_seq
python ./split_fastas.py Z_proteins.fa
mv *.fasta Z_seq ; cd Z_seq
for i in *.fasta ;  do blastp -query $i -db ../proteins_except_Zs -evalue 0.00001 -outfmt 6 -out $i_Z_blast_matches.txt; done
# once done:
rm list_of_input_files.txt

mkdir W_seq
python ./split_fastas.py W_proteins.fa 
mv *.fasta W_seq ;cd W_seq
ls *.fasta > list_of_input_files.txt
mbMem=50000; bsub -J "jobname[1-2000]%200" -o job_output.%J.%I.out -R"span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" -M${mbMem} ./run_blast_as_array.sh 
# for i in *.fasta ;  do blastp -query $i -db ../proteins_except_Ws -evalue 0.00001 -outfmt 6 -out $i_W_blast_matches.txt; done
# once done:
rm list_of_input_files.txt
cd ..

# run get_reciprocal_best_hits.py interactively
# results in a table that looks like this:
#anno1.g17206.t1 SUPER_Z1 anno1.g14927.t1 SUPER_W1
#anno1.g17671.t1 SUPER_Z2 anno1.g15496.t1 SUPER_W1
#anno1.g17236.t1 SUPER_Z1 anno1.g14883.t1 SUPER_W1
#anno1.g17343.t1 SUPER_Z2 anno1.g16250.t1 SUPER_W2
#anno1.g17663.t2 SUPER_Z2 anno1.g15486.t1 SUPER_W1
# saves list of reciprocal best hits to a tsv file called 'summarised_repicocal_blast_hits.txt'

grep SUPER_W1 summarised_repicocal_blast_hits.txt | grep SUPER_Z1 > summarised_repicocal_blast_hits_W1_Z1.txt # 156 hits
grep SUPER_W1 summarised_repicocal_blast_hits.txt | grep SUPER_Z2 > summarised_repicocal_blast_hits_W1_Z2.txt # 223 hits
grep SUPER_W2 summarised_repicocal_blast_hits.txt  | grep SUPER_Z1 > summarised_repicocal_blast_hits_W2_Z1.txt # 0 hits
grep SUPER_W2 summarised_repicocal_blast_hits.txt  | grep SUPER_Z2 > summarised_repicocal_blast_hits_W2_Z2.txt # 87 hits

# so going forward, only need to consider W1_Z1, W1_Z2, W2_Z2

# now lets extract their sequences and save to a folder ready for alignment
mkdir reciprocal_hits ; mkdir reciprocal_hits/proteins

# List of combinations
combinations=("W1_Z1" "W1_Z2" "W2_Z2")

# Second file path
second_file="BRAKER_1_BRAKER_2_TSEBRA.output.protein_sequences.oneliner.filtered.fa"

# Loop through each combination
for combination in "${combinations[@]}"; do
    input_file="summarised_repicocal_blast_hits_${combination}.txt"
    
    # Ensure the output directory exists
    output_dir="reciprocal_hits/proteins/${combination}"
    mkdir -p "$output_dir"
    
    # Read each line from the input file
    while IFS=' ' read -r col1 col2 col3 col4 col5; do
        output_file="${output_dir}/${col1}.plus_hit.fa"
        grep -A1 "$col1" "$second_file" > "$output_file"
        grep -A1 "$col3" "$second_file" >> "$output_file"
    done < "$input_file"
done
