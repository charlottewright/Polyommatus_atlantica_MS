
# Aim: Align each pair of protein sequences using FSA then translate the protein alignments into codon alignments pal2nal, and run with codeml

conda activate phylogenomics_env

cd reciprocal_hits/proteins/W1_Z1
for i in *.fa ; do fsa $i > `echo "$i" | cut -f1,2,3 -d'.'`.aligned.fna ; done
cd ../W1_Z2
for i in *.fa ; do fsa $i > `echo "$i" | cut -f1,2,3 -d'.'`.aligned.fna ; done
cd ../W2_Z2
for i in *.fa ; do fsa $i > `echo "$i" | cut -f1,2,3 -d'.'`.aligned.fna ; done
cd ../..
mkdir protein_alignments ; cd protein_alignments
mkdir W1_Z1 ; mkdir W1_Z2 ; mkdir W2_Z2 
mv ../proteins/W1_Z1/*.fna W1_Z1/
mv ../proteins/W1_Z2/*.fna W1_Z2/
mv ../proteins/W2_Z2/*.fna W2_Z2/

# download pal2nal - use to convert protein alignments into nucleotide alignments
cd ../../../../Scripts
curl -o pal2nal.pl https://raw.githubusercontent.com/LANL-Bioinformatics/PhaME/master/src/pal2nal.pl
chmod +x pal2nal.pl
cd -

mamba deactivate && mamba activate AGAT_env
agat_sp_extract_sequences.pl --gff BRAKER_1_BRAKER_2_TSEBRA.output.filtered.gtf -o BRAKER_1_BRAKER_2_TSEBRA.output.protein_nucleotide_sequences.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' BRAKER_1_BRAKER_2_TSEBRA.output.protein_nucleotide_sequences.fa > BRAKER_1_BRAKER_2_TSEBRA.output.protein_nucleotide_sequences.one_liner.fa

# split into one file per cdna sequence (needed as input to pal2nal)
python split_fastas.py BRAKER_1_BRAKER_2_TSEBRA.output.protein_nucleotide_sequences.one_liner.fa
mkdir nucleotide_sequences
mv anno* nucleotide_sequences/

# need to get specifically pairs of hits used in folders above
mkdir reciprocal_hits/ 
mkdir reciprocal_hits/nucleotide_sequences
mkdir reciprocal_hits/nucleotide_sequences/W1_Z1
mkdir reciprocal_hits/nucleotide_sequences/W1_Z2
mkdir reciprocal_hits/nucleotide_sequences/W2_Z2
#mv anno* nucleotide_sequences//W1_Z2 ;  mkdir reciprocal_hits/cdna_sequences/W2_Z2
for i in reciprocal_hits/protein_alignments/W1_Z1/* ; do ./get_corresponding_cdna_sequences.sh $i ; done
mv *plus_hit.fna reciprocal_hits/nucleotide_sequences/W1_Z1
for i in reciprocal_hits/protein_alignments/W1_Z2/* ; do ./get_corresponding_cdna_sequences.sh $i ; done
mv *plus_hit.fna reciprocal_hits/nucleotide_sequences/W1_Z2
for i in reciprocal_hits/protein_alignments/W2_Z2/* ; do ./get_corresponding_cdna_sequences.sh $i ; done
mv *plus_hit.fna reciprocal_hits/nucleotide_sequences/W2_Z2

# ready to do pal2nal1
mkdir reciprocal_hits/nucleotide_sequences_aligned
cd reciprocal_hits/nucleotide_sequences_aligned
mkdir W1_Z1 ; mkdir W1_Z2 ; mkdir W2_Z2

#pal2nal.pl reciprocal_hits/protein_alignments/W1_Z1/anno1.g16999.t1.aligned.fna reciprocal_hits/nucleotide_sequences/W1_Z1/anno1.g16999.t1.plus_hit.fna  -output paml > test.fa
# nb default output format is "clustal" which isnt recognised by codeml so specify either fasta or paml output when running pal2nal
for i in reciprocal_hits/protein_alignments/W1_Z1/*.fna ; do pal2nal.pl $i reciprocal_hits/nucleotide_sequences/W1_Z1/`echo "$i" | cut -f 4 -d'/' | cut -f1,2,3 -d'.'`.plus_hit.fna -output fasta > reciprocal_hits/nucleotide_sequences_aligned/W1_Z1/`echo "$i" | cut -f 4 -d'/' |  cut -f1,2,3,4 -d'.'`.aln ; done
for i in reciprocal_hits/protein_alignments/W1_Z2/*.fna ; do pal2nal.pl $i reciprocal_hits/nucleotide_sequences/W1_Z2/`echo "$i" | cut -f 4 -d'/' | cut -f1,2,3 -d'.'`.plus_hit.fna -output fasta > reciprocal_hits/nucleotide_sequences_aligned/W1_Z2/`echo "$i" | cut -f 4 -d'/' |  cut -f1,2,3,4 -d'.'`.aln ; done
for i in reciprocal_hits/protein_alignments/W2_Z2/*.fna ; do pal2nal.pl $i reciprocal_hits/nucleotide_sequences/W2_Z2/`echo "$i" | cut -f 4 -d'/' | cut -f1,2,3 -d'.'`.plus_hit.fna -output fasta > reciprocal_hits/nucleotide_sequences_aligned/W2_Z2/`echo "$i" | cut -f 4 -d'/' |  cut -f1,2,3,4 -d'.'`.aln ; done

# to run codeml, you first need to make a control file per alignment
# lets make a template ctrl file then edit it for each alignment
mkdir codeml ; cd codeml
touch test.nwk
mkdir W1_Z1
mkdir W1_Z2
mkdir W2_Z2
cp test.nwk W1_Z1/
cp test.nwk W1_Z2/
cp test.nwk W2_Z2/

# List of combinations
combinations=("W1_Z1")
# "W1_Z2" "W2_Z2")

# Loop through each combination
for combination in "${combinations[@]}"; do

    for i in ../reciprocal_hits/nucleotide_sequences_aligned/${combination}/*.aln; do 
        ctrl_file=`echo "$i" | cut -f5 -d'/' | cut -f1,2,3 -d'.'`.ctrl
        output_file=`echo "$i" | cut -f5 -d'/' | cut -f1,2,3 -d'.'`.codeml
        cp template_codeml.ctrl "$ctrl_file"
        sed -i "s|aln|$i|g" "$ctrl_file" # use '|' instead of '/' so sed doesn't confuse the '/' in the filepath with the sed command
        sed -i "s/output/$(basename "$output_file")/g" "$ctrl_file"
        codeml "$ctrl_file"
        mv *.ctrl ${combination}/
        mv *.codeml ${combination}/
        mv ${combination}/template_codeml.ctrl .
    done

done

cd W1_Z1
for i in *.ctrl ; do codeml $i; done

# summarise results
cd W1_Z1
for i in *.codeml ; do grep -B1 pairwise $i | grep anno >> ../W1_Z1.codeml_results.txt; done
for i in *.codeml ; do grep "t=" $i >> ../W1_Z1.codeml_ML_dS_results.txt; done
cd ../W1_Z2
for i in *.codeml ; do grep -B1 pairwise $i | grep anno >> ../W1_Z2.codeml_results.txt; done
for i in *.codeml ; do grep "t=" $i >> ../W1_Z2.codeml_ML_dS_results.txt; done
cd ../W2_Z2
for i in *.codeml ; do grep -B1 pairwise $i | grep anno >> ../W2_Z2.codeml_results.txt; done
for i in *.codeml ; do grep "t=" $i >> ../W2_Z2.codeml_ML_dS_results.txt; done
cd ..

# there are two inferred values of ds - one is from the NG method (quick, counting method)
awk '{print $4}' W1_Z1.codeml_NG_dS_results.txt | sed 's/)//' 
awk '{print $4}' W1_Z2.codeml_NG_dS_results.txt | sed 's/)//' 
awk '{print $4}' W2_Z2.codeml_NG_dS_results.txt | sed 's/)//' 

awk '{print $1, $4}' W1_Z1.codeml_NG_dS_results.txt | sed 's/)//' > W1_Z1.codeml_NG_dS_results.reformatted.txt
awk '{print $1, $4}' W1_Z2.codeml_NG_dS_results.txt | sed 's/)//' > W1_Z2.codeml_NG_dS_results.reformatted.txt
awk '{print $1, $4}' W2_Z2.codeml_NG_dS_results.txt | sed 's/)//' > W2_Z2.codeml_NG_dS_results.reformatted.txt

# the other ds value is from a max-likelihood approach (that takes multiple subs into account)
awk '{print $14}' W1_Z1.codeml_ML_dS_results.txt 
awk '{print $14}' W1_Z2.codeml_ML_dS_results.txt 
awk '{print $14}' W2_Z2.codeml_ML_dS_results.txt 

# get locations of transcripts
awk '$3 == "transcript"' BRAKER_1_BRAKER_2_TSEBRA.output.gtf | awk '{print $9, $1, $4, $5}' > BRAKER_1_BRAKER_2_TSEBRA.output.transcript_locations.tsv
