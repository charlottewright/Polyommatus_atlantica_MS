#!/usr/bin/env python3
#%%

#import matplotlib.pyplot as plt
import random
import argparse

def parse_P_atlantica_fai(fai_file):
    chrom_sizes, prop_chrom_sizes = [], []
    chrom_sizes_dict = {}
    total_genome_size = 0
    # Open the .fai file
    with open(fai_file, 'r') as file:
        for line in file:
            # Split each line by tab
            columns = line.strip().split('\t')
            chrom_name = columns[0]
            chrom_size = int(columns[1])
            
            # Add to total genome size
            if 'SUPER_W' not in chrom_name: # to get haploid genome_size (i.e don't include Ws as well as Zs)
                total_genome_size += chrom_size

            # if an autosome, add to list of chrom_sizes
            if ('SUPER' in chrom_name) and ('unloc' not in chrom_name) and ('SUPER_W' not in chrom_name) and ('SUPER_Z' not in chrom_name):
                chrom_sizes.append(chrom_size)
                chrom_sizes_dict[chrom_name] = chrom_size
    for i in chrom_sizes:
        prop_chrom_sizes.append(i/total_genome_size)
    return total_genome_size, chrom_sizes, chrom_sizes_dict, prop_chrom_sizes

def parse_P_icarus_fai(fai_file):
    chrom_sizes, prop_chrom_sizes = [], []
    chrom_sizes_dict = {}
    total_genome_size = 0
    # Open the .fai file
    with open(fai_file, 'r') as file:
        for line in file:
            # Split each line by tab
            columns = line.strip().split('\t')
            chrom_name = columns[0]
            chrom_size = int(columns[1])
            
            # Add to total genome size
            total_genome_size += chrom_size

            # if an autosome, add to list of chrom_sizes (i.e. not Z or mito)
            if ('OW56' in chrom_name) and ('OW569320.1' not in chrom_name) and ('OW569343.1' not in chrom_name):
                chrom_sizes.append(chrom_size)
                chrom_sizes_dict[chrom_name] = chrom_size
        for i in chrom_sizes:
            prop_chrom_sizes.append(i/total_genome_size)
    return total_genome_size, chrom_sizes, chrom_sizes_dict, prop_chrom_sizes

def break_chromosomes_old(chromosomes, target_lengths, min_length=2, max_retries=1):
    result = {}    # dictionary to hold the resulting chromosomes pieces after breaking
    for chr, length in chromosomes.items():
        if len(target_lengths) == 0:
            print('Cannot fragment chr as we are out of lengths!')
        current_length = length
        current_chr = chr
        pieces = []  # To store the broken pieces of the current chr
        remaining_target_lengths = target_lengths#.copy()  # Copy the target lengths for each chr
        start_position = 1  # Start at position 1 for each chr
        
        while current_length >= min_length and remaining_target_lengths:
            random.shuffle(remaining_target_lengths)  # Shuffle to randomly pick from target_lengths
            retries = 0
            piece_found = False
            
            # Try breaking the chr based on target lengths, with retries in case of invalid break
            while retries < max_retries and remaining_target_lengths:
                target = remaining_target_lengths.pop()  # Pick a random target length
                target = int(target)
                # If the target is smaller than the remaining length and above the min_length, break
                #if target <= current_length and target >= min_length:
                if target <= current_length and (current_length-target) >= min_length:
                    pieces.append((current_chr, start_position, start_position + target -1))
                    start_position += target  # Update the start position for the next piece
                    current_length -= target  # Reduce the remaining length of the chr
                    current_chr = f'{chr}_{len(pieces) + 1}'  # Create a new chr name for the next part
                    piece_found = True
                    break
                else: # if breaking a chr was unsuccessful, add this length back into the pack for use by other chr (I think)
                    remaining_target_lengths.append(target)
                retries += 1
            
            # If no valid break was found after max_retries, stop trying to break this chr
            if not piece_found:
                print(f"Failed to break {chr} further after {max_retries} retries.")
                break

        # After all valid breaks, if there is any remaining length, make it the last piece
        if current_length > 0:
            pieces.append((current_chr, start_position, start_position + (current_length-1)))

        # Store the resulting pieces in the result dictionary
        result[chr] = pieces

    return result

def write_bed_file(tuples_list, filename='output.bed', outputdir='./'):
    counter = 1
    if filename != 'simulated_chr.bed':
            filename = filename + '.simulated_chr.bed'
    filepath = outputdir + filename
    with open(filepath, 'w') as f:
        print('File written to:', filepath)
        for piece in tuples_list:
            chr = piece[0].split('_')[0]
            sim_chr_ID = "SimChr" + str(counter)
#            f.write(f"{chr} {piece[1]} {piece[2]}\n")
            f.write(f"{chr}\t{piece[1]}\t{piece[2]}\t{sim_chr_ID}\n")
            counter = counter + 1


def break_chromosomes(P_icarus_chrom_sizes_dict, target_chr_lengths, min_length):
    P_icarus_starting_chrom_sizes_dict = {} # keeps track of starting points for each new fragment
    for i in P_icarus_chrom_sizes_dict.keys():
        P_icarus_starting_chrom_sizes_dict[i] = 0
    P_icarus_starting_chrom_sizes_dict

    chrom_count = {}
    pieces = []
    number_fragments_generated = 0
    target_chr_lengths_copy = target_chr_lengths.copy()
    P_icarus_chrom_sizes_dict_copy = P_icarus_chrom_sizes_dict.copy() # keeps track of remaining chr length
    P_icarus_chr_list = list(P_icarus_chrom_sizes_dict_copy.keys())

    # if the random_length is less than the remaining chr length
    # and the remaining length after the random length results in a fragment at least the size of the min length:
    while len(pieces) < 205:
        random.shuffle(target_chr_lengths_copy)  # Shuffle to randomly pick from target_lengths
        random.shuffle(P_icarus_chr_list)  # Shuffle to randomly pick from chr
        random_length = target_chr_lengths_copy.pop()  # Pick a random target length
        random_chr = P_icarus_chr_list.pop()  # Pick a random target length
        current_chr_length = P_icarus_chrom_sizes_dict_copy[random_chr]

        if random_length <= current_chr_length and (current_chr_length-random_length >= min_length):
                            # lets find a number to label this piece
                            if random_chr in chrom_count.keys():
                                piece_number = chrom_count[random_chr] + 1
                                chrom_count[random_chr] = piece_number # update value
                            else: 
                                piece_number = 1
                                chrom_count[random_chr] = 1
                            new_piece = f'{random_chr}_{piece_number}'  # Create a new chr name for the next part
                            # look up start position
                            start_pos = P_icarus_starting_chrom_sizes_dict[random_chr]
                            pieces.append((new_piece, start_pos, start_pos+random_length))
                            # now updat the starting position in the dict
                            P_icarus_starting_chrom_sizes_dict[random_chr] = start_pos+random_length
                            # and update remaining chromosome length
                            P_icarus_chrom_sizes_dict_copy[random_chr] = current_chr_length-random_length
                            number_fragments_generated = number_fragments_generated+1
        else:
         #   print('Fragment was not permitted as it either it would have created too small fragments or the chromosome length was less than the fragment size')
            target_chr_lengths_copy.append(random_length) # add random length back into the list
        P_icarus_chr_list.append(random_chr) # add chr back to chr list

    # now lets assign the residual pieces of each chr (i.e. remaining lengths)
    for chr, starting_length in P_icarus_starting_chrom_sizes_dict.items():
        remaining_length = P_icarus_chrom_sizes_dict_copy[chr]
        chr_end = starting_length+remaining_length
        if chr in chrom_count.keys():
            piece_number = chrom_count[chr] + 1
            chrom_count[chr] = piece_number # update value
        else: 
            piece_number = 1
            chrom_count[chr] = 1
        new_piece = f'{chr}_{piece_number}'  # Create a new chr name for the next part
        # look up start position
        pieces.append((new_piece, starting_length, chr_end))

    return(pieces, chrom_count)
#%%

# Initialize the parser
parser = argparse.ArgumentParser(description="Simulate_fissions_by_length.py")

# Define the three string arguments
parser.add_argument('--ref_index', type=str, required=True, help="Index file of species which has fissions")
parser.add_argument('--query_index', type=str, required=True, help="Index file of species to be used to simulate fissions")
parser.add_argument('--prefix', type=str, required=True, help="output prefix for output file")
parser.add_argument('--output_dir', type=str, help="output directory for output file", default='./')

# Parse the command-line arguments
args = parser.parse_args()

# Print out the arguments
print(f"Reference index (genome with fissions): {args.ref_index}")
print(f"Query index (genome used to simulate fissions): {args.query_index}")
print(f"Directory for output file: {args.output_dir}")
print(f"Prefix for output bed file: {args.prefix}")
#%%
P_atlantica_fai_file = args.ref_index
P_icarus_fai_file = args.query_index
output_prefix = args.prefix
output_dir = args.output_dir
#P_atlantica_fai_file = '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Raw_data/genomes/core_dataset/Polyommatus_atlantica.fai'  # Path to .fai file
#P_icarus_fai_file = '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Raw_data/genomes/core_dataset/Polyommatus_icarus.fa.fai'  # Path to .fai file

P_atlantica_genome_size, P_atlantica_chromosome_sizes, P_atlantica_chr_sizes_dict, P_atlantica_proportional_chromosome_sizes = parse_P_atlantica_fai(P_atlantica_fai_file)

# Print results
print(f"Genome size of reference genome (genome with fissions): {round(P_atlantica_genome_size/1000000, 2)} Mb")
print(f"Number of chromosomes in reference genome: {len(P_atlantica_chromosome_sizes)}")
print(f"Minimum proportional length of a chromosome in the reference genome: {round(min(P_atlantica_proportional_chromosome_sizes),2)}")
# %%
P_icarus_genome_size, P_icarus_chromosome_sizes, P_icarus_chrom_sizes_dict, P_icarus_proportional_chromosome_sizes = parse_P_icarus_fai(P_icarus_fai_file)

# set minimum length to be the smallest chr in P. atlantica scaled by genome size in P. icarus
min_length = int(min(P_atlantica_proportional_chromosome_sizes)*P_icarus_genome_size)
target_chr_lengths = [int(i*P_icarus_genome_size) for i in P_atlantica_proportional_chromosome_sizes]

#%%
P_atlantica_proportional_chromosome_sizes
#%%
pieces, chrom_count = break_chromosomes(P_icarus_chrom_sizes_dict, target_chr_lengths, min_length)

# %%
# lets plot the distribution of fragments per P. icarus chr
# Known example: chr 7 (OW569327.1) is lenth 24 Mb and splits into 14 chr 
# Largest chr: chr 1 (OW569321.1) should split into ~19
# Smallest chr: chr 22 (OW569342.1) should split into 5 fragments
# Sorting the keys of both dictionaries
sorted_keys = sorted(P_icarus_chrom_sizes_dict.keys())  
sorted_counts = [chrom_count[key] for key in sorted_keys]
sorted_lengths = [P_icarus_chrom_sizes_dict[key] for key in sorted_keys]

# save results
sorted_pieces = sorted(pieces, key=lambda x: (x[0], int(x[1])))
write_bed_file(sorted_pieces, output_prefix, output_dir)

#write_bed_file(sorted_pieces, '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus_atlantica/Analysis/enrichment_analysis/P_atlantica_enrichment/single_copy_orthos_with_ilPolIcar1/test2.bed')

# %%
# Create a dot plot
#plt.scatter(sorted_lengths, sorted_counts, color='blue')

# Add labels and title
#plt.xlabel('P. icarus chr length (Mb)')
#plt.ylabel('Number of fragments')
#plt.title('Result of breaking chr')

# Show the plot
#plt.show()

#%%
