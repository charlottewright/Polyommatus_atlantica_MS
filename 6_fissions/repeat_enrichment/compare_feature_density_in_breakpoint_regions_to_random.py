
#!/usr/bin/python3

# Aim: ask whether repeat density is higher/lower than in rest of chr on average


#breakpoint_regions<-import(args[1])
#repeats <-import(args[2])


# numOverlaps(breakpoint_regions, repeats, count.once=TRUE)

# number of times per breakpoint that overlaps with repeats
# ask is this higher or lower tha
# overlapPermTest(A=breakpoint_regions, B=repeats, ntimes=as.numeric(args[6]), genome=genome, non.overlapping=FALSE, count.once=TRUE)
#%%
import sys
import argparse
import statistics
import seaborn as sns
import matplotlib.pyplot as plt
import math
import random

#%%
def read_breakpoints_file(breakpoints_file):
    chr2breakpoint = []
    with open(breakpoints_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            chr, start, end, density = cols[0], int(cols[3]), int(cols[4]), float(cols[6])
            breakpoint = (chr, start, end, density)
            chr2breakpoint.append(breakpoint)
    return(chr2breakpoint)

def read_random_regions_file(file):
    chr2repeat_density = {}
    with open(file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            chr, density = cols[0], float(cols[6])
            if chr in chr2repeat_density.keys():
                current_density_values = chr2repeat_density[chr]
                current_density_values.append(density)
                
            else:
                chr2repeat_density[chr] = [density]
    return(chr2repeat_density)

# implement function used in regioneR https://www.biostars.org/p/431716/ to calculate p-value and z-score
def get_p_value(random_values, observed_value, side):
    if side == "greater": # alternative hypothesis is that the observed value is greater than null distribution
        filtered_values = [i for i in random_values if i >= observed_value] # if random value is greater than observed value
    else:
        filtered_values = [i for i in random_values if i <= observed_value] # if random value is less than observed value
    pval = (len(filtered_values) + 1) / (len(random_repeat_density_values) + 1)
    return(pval)

def draw_histogram(random_values, observed_value): # e.g. random_repeat_density_values, breakpoint_density
    plt.figure()
    plot_title = 'Repeat density across random regions of chromosome:' + chr
    hist_plot = sns.histplot(data=random_values)
    hist_plot.axvline(x = observed_value, color="red")
    hist_plot.set(title=plot_title)
    return()

def write_output(pval_less_list, pval_greater_list, total_breakpoints):
    average_less_than_pval = sum(pval_less_list)/total_breakpoints
    average_greater_than_pval = sum(pval_greater_list)/total_breakpoints
    average_zscore = sum(zscore_list)/total_breakpoints
    with open(output_file, 'w') as file:
        file.write("%s\t%s\t%s\t%s" % ('comparison', 'average_less_than_pval', 'average_greater_than_pval', 'average_zscore') + "\n")
        file.write("%s\t%s\t%s\t%s" % (output_prefix, average_less_than_pval, average_greater_than_pval, average_zscore) + "\n")
    return()
#%%

if __name__ == "__main__":
    SCRIPT = "random_region_generator.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--breakpoints_file", type=str, help = "locations of breakpoint regions and their repeat densities in TSV format", required=True)
    parser.add_argument("-o", "--output_prefix", type=str, help = "prefix for output file names", default="")
    parser.add_argument("-r", "--random_regions_file", type=str, help = "locations of random regions and their repeat densities in TSV format", required=True)
    args = parser.parse_args()
    breakpoints_file = args.breakpoints_file
    output_prefix = args.output_prefix
    random_regions_file = args.random_regions_file
    
    chr2breakpoint = read_breakpoints_file(breakpoints_file)
    chr2random_repeat_density = read_random_regions_file(random_regions_file)

#    breakpoints_file = 'repeat_enrichment/repeat_density_of_ilCyaSemi1.1_vs_ilLysCori1.1_breakpoints_refined_231123.tsv_Cyaniris_semiargus.filteredSINE.gff.tsv'
#    random_regions_file = 'repeat_enrichment/repeat_density_of_random_regions_10000_Cyaniris_semiargus.filteredSINE.gff.tsv'

    # only thing I'm not sure about is how to take into account different chr - i am calculating z-score & p-value per breakpoint, then taking average z-score and p-value - might be wrong?
    print('[+]      Calculating p-value and Z-score per breakpoints:')

    pval_greater_list, pval_less_list = [], []
    zscore_list = []
    total_breakpoints = len(chr2breakpoint)
    for breakpoint_info in chr2breakpoint:
        num_greater = 0 # number of times the random is greater than observed
        num_less = 0
        num_zeros = 0
        chr = breakpoint_info[0]
        breakpoint_density = breakpoint_info[3]
        random_repeat_density_values = chr2random_repeat_density[chr]
        less_pval = get_p_value(random_repeat_density_values, breakpoint_density, "less")
        greater_pval = get_p_value(random_repeat_density_values, breakpoint_density, "greater")
        zscore = round((breakpoint_density - statistics.mean(random_repeat_density_values))/statistics.stdev(random_repeat_density_values), 4)
        pval_less_list.append(less_pval)
        pval_greater_list.append(greater_pval)
        zscore_list.append(zscore)
      #  print('Breakpoint density is:', breakpoint_density)
      #  print('Average random density is:', statistics.mean(random_repeat_density_values))
      #  print('Number of random points with density greater than breakpoint:', num_greater)
      #  print('Number of random points with density less or equal to the breakpoint:', num_less)
      #  print('Number of random regions with a density of zero:', num_zeros)
     #   print('Chr:', chr, ' Greater pval: ', round(greater_pval,2), ' Less pval: ', round(less_pval,2), ' Zscore: ', zscore)
    print('[+]      Generating average p-values and Z-scores across all breakpoints:')
    print('Average less than pval: ', sum(pval_less_list)/total_breakpoints)
    print('Average greater than pval: ', sum(pval_greater_list)/total_breakpoints)
    print('Average zscore: ', sum(zscore_list)/total_breakpoints)
    print('')

    for percentile in [0.01, 0.1, 1, 5, 10, 50, 90, 95, 99, 99.8]:
        repeat_enriched_breakpoints = 0
        for breakpoint_info in chr2breakpoint:
            chr = breakpoint_info[0]
            breakpoint_density = breakpoint_info[3]
            random_repeat_density_values_sorted = sorted(chr2random_repeat_density[chr])
            random_upper_bound = random_repeat_density_values_sorted[math.ceil(len(random_repeat_density_values) * (percentile / 100))]
            if breakpoint_density > random_upper_bound:
                repeat_enriched_breakpoints +=1
        # Summarise results
        print("{:.2f}".format((repeat_enriched_breakpoints/total_breakpoints)*100), '% of breakpoints have repeat density greater than the', percentile, '% percentile amongst random regions.')

output_file = 'breakpoint_vs_random_regions_' + output_prefix + '.tsv'
print(' [+]     Writing results to output ', output_file)
write_output(pval_less_list, pval_greater_list, total_breakpoints)

#%%

exit()
# Plot all regions together - NB this no longer shows the nuance of different repeat_densities per chr
random_repeat_density_values_total = []
breakpoint_densities = []
for breakpoint_info in chr2breakpoint:
        num_greater = 0 # number of times the random is greater than observed
        num_less = 0
        num_zeros = 0
        chr = breakpoint_info[0]
        breakpoint_density = breakpoint_info[3]
        random_repeat_density_values = chr2random_repeat_density[chr]
        random_repeat_density_values_total.extend(random_repeat_density_values)
        breakpoint_densities.append(breakpoint_density)
#%%
plt.figure()
plot_title = 'Repeat density across random regions vs in breakpoints on average'
random_repeat_density_values_subset = random.sample(random_repeat_density_values_total, k=1000)
hist_plot = sns.histplot(data=random_repeat_density_values_subset, bins=100)
hist_plot.axvline(x = statistics.mean(breakpoint_densities), color="red") # average repeat density in breakpoints
hist_plot.set(title=plot_title)
# %%
import pandas as pd
df = pd.DataFrame({'random_regions': random_repeat_density_values_total, 'breakpoints': breakpoint_densities})

# %%
# Compare repeat density distibutions between random vs breakpoints - more reflective of the different landscapes of repeat densities per chr
random_repeat_density_values_subset = random.sample(random_repeat_density_values_total, k=len(breakpoint_densities))
plt.hist(random_repeat_density_values_subset, bins=50, color='blue', alpha=0.7, label='List 1')
plt.hist(breakpoint_densities, bins=50, color='orange', alpha=0.7, label='List 2')
# %%
plt.figure()
plot_title = 'Distribution of pvals'
hist_plot = sns.histplot(data=pval_greater_list)
hist_plot.axvline(x = 0.05, color="red")
hist_plot.set(title=plot_title)
# %%
# multiple testing correction
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

# Create a list of the adjusted p-values
np.array(pval_greater_list) * len(pval_greater_list)
# or divide alpha value by len(pval_list) to get new cutoff for significance
adjusted_pvalue_cutoff = 0.05/total_breakpoints
significant_pvals = [i for i in pval_greater_list if i < adjusted_pvalue_cutoff]
print('Number of significant results:', len(significant_pvals))
