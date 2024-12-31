
import random

#%%
def get_locations_of_query_genes(gene_locations_file):
    chr2locations = {}
    with open(gene_locations_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                cols = line.split('\t')
                if cols[2] == "gene":
                    chr = cols[0]
                    start, stop = int(cols[3]), int(cols[4])
                    if start > stop:
                        temp = start
                        start = stop
                        stop = temp
                    if ('Z' not in chr) & ('CALM' not in chr): # i.e. exclude the Z chr and unplaced scaffs
                        if chr in chr2locations:
                            chr2locations[chr].append(stop)
                        else:
                            chr2locations[chr] = [stop]
    return(chr2locations)

def get_gene_counts (gene_counts_file):
    chr2count = {}
    with open(gene_counts_file, 'r') as file:
        for line in file:
            cols = line.split('\t')
            chr = cols[0]
            count = int(cols[3])
            length = int(cols[2])
            if (length >= 1200000) & ('W' not in chr) & ('Z' not in chr):
                chr2count[chr] = count
    gene_counts = [i for i in chr2count.values()]
    return(gene_counts)

def generate_new_fragments_old(chr2locations, gene_counts):
    new_chr = []   # list of tibbles of chr_ID, start, en

    for chr, locations in chr2locations.items():
        current_gene = 1  # start on the first location
        position_to_start_at = 1  # start first chr at start of a chr
        locations_sorted = sorted(locations)
        number_locations = len(locations)
        position_to_lookup = 1  # initialise at '1' at start of each chr loop

        while position_to_lookup <= number_locations:  # i.e., until we've run out of locations
            # Generate a random number of genes in a chr
            number_genes = random.choice(gene_counts)
        #   print('Random number of genes to add to pos:', number_genes)
            gene_counts.remove(number_genes)  # remove value from list so we don't draw it again
            position_to_lookup = current_gene + number_genes  # adjust lookup position
            if position_to_lookup >= number_locations:
                new_chr.append((chr, position_to_start_at, max(locations), number_genes)) # we're done with this chr so add an entry for rest of chr, as remember this will be a fragment too
                # If position_to_lookup exceeds available locations, exit the loop
                break
            position_to_end_at = locations_sorted[position_to_lookup]
        #   print('So the end of the new chr is:', position_to_end_at)
            # Assuming 'chr' is defined elsewhere or you may need to define it
        #  print('New chr:', chr, position_to_start_at, position_to_end_at)            
            new_chr.append((chr, position_to_start_at, position_to_end_at, number_genes))
            position_to_start_at = position_to_end_at  # set pos to start at the end of the last cut
            current_gene = position_to_lookup + 1  # move to the next gene position
    return(new_chr, gene_counts)

def generate_new_fragments(chr2locations, scaled_gene_count_dist):
    new_chr = []   # list of tibbles of chr_ID, start, en
    genes_per_new_chr = [] # lets keep track of genes per new chr
    for chr, locations in chr2locations.items():
        current_gene = 1  # start on the first location
        position_to_start_at = 1  # start first chr at start of a chr
        locations_sorted = sorted(locations)
        number_locations = len(locations)
        position_to_lookup = 1  # initialise at '1' at start of each chr loop

        while position_to_lookup < number_locations:  # i.e., until we've run out of locations
            # Generate a random number of genes in a chr
            number_genes = random.choice(scaled_gene_count_dist)
        #   print('Random number of genes to add to pos:', number_genes)
            position_to_lookup = current_gene + number_genes  # adjust lookup position
            if position_to_lookup >= number_locations:
                new_chr.append((chr, position_to_start_at, max(locations), number_genes)) # we're done with this chr so add an entry for rest of chr, as remember this will be a fragment too
                number_genes_in_new_chr = len(locations) - current_gene
                genes_per_new_chr.append(number_genes_in_new_chr)
                # If position_to_lookup exceeds available locations, exit the loop
                break
            else:
                scaled_gene_count_dist.remove(number_genes)  # remove value from list so we don't draw it again
                position_to_end_at = locations_sorted[position_to_lookup]
        #   print('So the end of the new chr is:', position_to_end_at)
            # Assuming 'chr' is defined elsewhere or you may need to define it
        #  print('New chr:', chr, position_to_start_at, position_to_end_at)            
                new_chr.append((chr, position_to_start_at, position_to_end_at, number_genes))
                genes_per_new_chr.append(number_genes)
                position_to_start_at = position_to_end_at  # set pos to start at the end of the last cut
                current_gene = position_to_lookup + 1  # move to the next gene position
    return(new_chr, scaled_gene_count_dist,genes_per_new_chr)

#%%
# get set of gene counts for P. atlantica per autosome
gene_counts_file = '../Analysis/gene_prediction/analysis/Polyommatus_atlantica.braker1_and_braker2.TSEBRA.gene_density_per_scaff.tsv'
gene_count_dist = get_gene_counts(gene_counts_file)

# read in locations of P. icarus genes
gene_locations_file = '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Analysis/enrichment_analysis/Polyommatus_icarus-GCA_937595015.1-2022_06-genes.filtered.modified.gff3'
chr2locations = get_locations_of_query_genes(gene_locations_file)

#%%
print('Total number of genes in P. atlantica:', sum(gene_count_dist))  # 16542
print('Average number of genes per autosome in P. atlantica:', round(sum(gene_count_dist)/len(gene_count_dist)))
print('Number of autosomes in P. atlantica:', len(gene_count_dist))

# started with 23 autosomes, ended with 227, so number of fission events was 204
# so need 204 entries in list to pick from
#print('number of values to remove from list:', 227-204)
print('Final starting set of gene counts has length of:', len(gene_count_dist))

total_genes = []
for chr, locs in chr2locations.items():
    for i in locs:
        total_genes.append(i)
# total num genes on autosomes in P. icarus
print('Number of genes in P. icarus:', len(total_genes)) # 13345
print('Number of chr in P. icarus:', len(chr2locations.keys())) # 13345
# %%
# total num genes on autosomes in P. atlantica
gene_count_dist = get_gene_counts(gene_counts_file)
num_genes_P_atlantica = sum(gene_count_dist)
num_genes_P_icarus = len(total_genes) # 13345
prop_difference = num_genes_P_icarus/num_genes_P_atlantica # P. icarus annotation has 25% less genes that P. atlantica
scaled_gene_count_dist = [] # lets adjust each number to be be consistent with numbers of genes in P. icarus
for i in gene_count_dist:
    scaled_value = int(i*prop_difference)
    scaled_gene_count_dist.append(scaled_value)

# %%
import matplotlib.pyplot as plt
import numpy as np

plt.hist(gene_count_dist)
plt.show() 

plt.hist(scaled_gene_count_dist)
plt.show() 

# %%
# First lets see the distribution of chr numbers we get when using this set
count = 0
simulated_chr_sets = []
while (count < 50):
#while (count < 1001):
    gene_count_dist = get_gene_counts(gene_counts_file)
    for i in gene_count_dist:
        scaled_value = int(i*prop_difference)
        scaled_gene_count_dist.append(scaled_value)
  #  for i in range(1, 24):
   #     val_remove = random.choice(scaled_gene_count_dist) # pick a value at random
    #    scaled_gene_count_dist.remove(val_remove) # remove value from list
    new_chr_set, final_counts = generate_new_fragments(chr2locations, scaled_gene_count_dist)
    simulated_chr_sets.append(new_chr_set)
    count = count + 1

#%%
len_dist = []
for i in simulated_chr_sets:
    len_dist.append(len(i))
plt.hist(len_dist)
plt.show() 
# conclude looks reasonably normally distributed, lets go ahead and restrict to 227 chromosomes
#%%
# now lets actually run the simulation - this time resticting results to those with 227 chromosomes
count = 0
simulated_chr_sets = []
while (count < 1000): # lets run 1000 simulations
    gene_count_dist = get_gene_counts(gene_counts_file)
    for i in gene_count_dist:
        scaled_value = int(i*prop_difference)
        scaled_gene_count_dist.append(scaled_value)
  #  for i in range(1, 24):
   #     val_remove = random.choice(scaled_gene_count_dist) # pick a value at random
    #    scaled_gene_count_dist.remove(val_remove) # remove value from list
    new_chr_set, final_counts, genes_per_new_chr = generate_new_fragments(chr2locations, scaled_gene_count_dist)
    if len(new_chr_set) == 227:
        simulated_chr_sets.append(new_chr_set)
        count = count + 1

#%%
plt.hist(genes_per_new_chr)
# looks similar to distribution of scaled genes per chr i.e. simulation is working as expected
#%%
len(simulated_chr_sets)
#%%
chr_set_ID = 1
output_simulations_file = '../Analysis/enrichment_analysis/P_icarus_similated_sets_227_autosomes.0824.tsv'

with open(output_simulations_file, 'w') as file:
    file.write("%s\t%s\t%s\t%s\t%s" % ("chr_set_ID", "chr_set", "start", "end", "number_genes") + "\n")
    for chr_set in simulated_chr_sets:
        for i in chr_set:
            file.write("%s\t%s\t%s\t%s\t%s" % (chr_set_ID, i[0], i[1], i[2], i[3]) + "\n")
        chr_set_ID = chr_set_ID + 1
 
