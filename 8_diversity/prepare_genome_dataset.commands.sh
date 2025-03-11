#!/usr/bin/bash

# NB this will download all Lycaenidae genomes that are publically available as of 09/07/24 so will not include ilPolAtla - this will be added to the dir separately.
cd /lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Raw_data/genomes/metadata/
mamba activate ncbi_datasets
datasets download genome taxon "lycaenidae" --reference --dehydrated --filename lycaenidae_090724.zip
dataformat tsv genome --package lycaenidae_090724.zip --fields accession,organism-name,assmstats-scaffold-n50,assminfo-level > all_lycaenidae_genomes_090724.tsv
sed -i 's/ /_/g' all_lycaenidae_genomes_090724.tsv 
mamba deactivate && mamba activate Basic_python
python3 ./filter_records.py all_lycaenidae_genomes_090724.tsv filtered_lepidoptera_genomes_09072024.tsv

accession="$(cut -f1 filtered_lepidoptera_genomes_09072024.tsv | tail +2| tr '\n' ' ')"
mamba deactivate && mamba activate ncbi_datasets
datasets download genome accession $accession --dehydrated --filename lycaenidaeFiltered.zip

unzip lycaenidaeFiltered.zip 
mkdir ncbi_dataset/ncbi_dataset
cp ncbi_dataset/fetch.txt  ncbi_dataset/ncbi_dataset/
datasets rehydrate --directory ncbi_dataset # this will take a little while

mv ncbi_dataset/ncbi_dataset/data/* ../
tail +2 filtered_lepidoptera_genomes_09072024.tsv | cut -f1,2 > ../GCA_accession_ID_2_species_name.tsv
cd ..
mv GC*/*.fna .
mv *.fa core_dataset/
# First move the male version of P. icarus to 
while IFS=$'\t' read -r GCA species_name ; do cp "${GCA}"*.fna "${species_name}.fa" ; done < GCA_accession_ID_2_species_name.tsv

# clean up
mv *fna genome_copies/
mv GCA_accession_ID_2_species_name.tsv metadata/
rm -R GC*
tar czf genome_copies.tar.gz genome_copies/
