#!/bin/bash

topgo_dir="$1"

# enter the topgo results dir
cd "$topgo_dir"

# Loop through each .txt file in the directory
for i in Sim*.BP; do
	awk '{print $3}' $i | sed 's/\"//g' | sed -n '2p' >> BP_terms_per_fragment.tsv
done


# save the mean number of BP terms per fragment for this simulated set to a summary file
awk 'BEGIN{s=0;}{s+=$1;}END{print s/NR;}' BP_terms_per_fragment.tsv >> ../mean_BP_per_fragment_per_simulation_set.03032025.tsv

# save the standard deviation of BP terms per fragment for this simulated set to a summary file
awk '{sum+=$0;a[NR]=$0}END{for(i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}' BP_terms_per_fragment.tsv  >> ../stdev_BP_per_fragment_per_simulation_set.03032025.tsv

# go back to starting directory
cd ..

