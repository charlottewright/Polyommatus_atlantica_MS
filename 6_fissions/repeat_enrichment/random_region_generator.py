
#!/usr/bin/env python3
import sys
import argparse
# Aim: ask whether repeat density is higher/lower than in rest of chr on average


#breakpoint_regions<-import(args[1])
#repeats <-import(args[2])


# numOverlaps(breakpoint_regions, repeats, count.once=TRUE)

# number of times per breakpoint that overlaps with repeats
# ask is this higher or lower tha
# overlapPermTest(A=breakpoint_regions, B=repeats, ntimes=as.numeric(args[6]), genome=genome, non.overlapping=FALSE, count.once=TRUE)
#%%
import random
from selectors import EpollSelector
import statistics
import seaborn as sns
import matplotlib.pyplot as plt

#%%


def read_breakpoints_file(breakpoints_file):
    chr2breakpoint = []
    with open(breakpoints_file, 'r') as file:
        for line in file:
            if 'start' not in line: # skip header
                cols = line.strip().split('\t')
                chr, start, end = cols[0], int(cols[1]), int(cols[2]) # edited
                breakpoint = (chr, start, end)
                chr2breakpoint.append(breakpoint)
    return(chr2breakpoint)

def calculate_mean_breakpoint_size(chr2breakpoint):
    breakpoint_sizes = []
    for i in chr2breakpoint:
        breakpoint_size = i[2] - i[1]
        breakpoint_sizes.append(breakpoint_size)
    return(statistics.mean(breakpoint_sizes), statistics.stdev(breakpoint_sizes))

def generate_random_region(mean_difference, std_deviation, chr_length):
    # Generate random start number greater than min_start
    max_start = int(chr_length - mean_difference - std_deviation)
    start_number = random.randint(0, max_start)
    while (start_number < 0):
        start_number = random.randint(0, max_start)
    # Generate random end number less than max_end
    end_number = round(start_number + abs(random.normalvariate(mean_difference, std_deviation)))
    while end_number > chr_length:
        end_number = round(start_number + abs(random.normalvariate(mean_difference, std_deviation)))
    result = (start_number, end_number)
    return result

def read_index(index_file):
    chr2length = {}
    with open(index_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            chr, length= cols[0], int(cols[1])
            if length < 50000: # use a filter of 50kb to remove mito and unplaced scaffolds
                    print('WARNING: contig ', chr, 'was removed due to its size being < 50 kb in length.')
            else:
                chr2length[chr] = length
    return(chr2length)

#%%
## Script ##
if __name__ == "__main__":
    SCRIPT = "random_region_generator.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--breakpoints_file", type=str, help = "locations of breakpoit regions in BED format", required=True)
    parser.add_argument("-i", "--index_file", type=str, help = "fasta index file of genome.", required=True)
    parser.add_argument("-o", "--output_prefix", type=str, help = "prefix for output file names", default="")
    parser.add_argument("-n", "--number_iterations", type=int, help = "number of random regions per chromosome to generate", default=1000)
    args = parser.parse_args()
    index_file = args.index_file
    breakpoints_file = args.breakpoints_file
    output_prefix = args.output_prefix
    number_iterations = args.number_iterations
#    breakpoints_file = '
#    index_file = 'Cyaniris_semiargus.fasta.fai'
#    total_sims = 100
#    output_file = 'test_output.tsv'

    # Read in fasta index file
    chr2length = read_index(index_file)
    # Read in breakpoints
    chr2breakpoint = read_breakpoints_file(breakpoints_file)
    # get average/stdev breakpoint size (across all breakpoints) to use to simulate random regions
    mean_breakpoint_size, stdev_breakpoint_size = calculate_mean_breakpoint_size(chr2breakpoint)
    # set these values here
    mean_difference = mean_breakpoint_size # 20000
    std_deviation = stdev_breakpoint_size # 10

    # make list of all chr with at least one breakpoint
    chr_list = list(set([item[0] for item in chr2breakpoint]))

    # Generate random regions of length of approx size of breakpoints. Make a set per chr.
    output_file = 'locations_of_random_regions_' + str(number_iterations) + '_per_chr_' + output_prefix + '.tsv'
    print(' [+]     Generating ', number_iterations, ' random regions per chr which are the approx size of breakpoints.' )
    with open(output_file, 'w') as file:
        for chr, chr_length in chr2length.items():
            for i in range(1,number_iterations+1, 1):
                random_region = generate_random_region(mean_difference, std_deviation, chr_length)
                file.write("%s\t%s\t%s" % (chr, random_region[0], random_region[1]) + "\n")
    print(' [+]     Locations have been written the the output file', output_file, 'in BED format.')
