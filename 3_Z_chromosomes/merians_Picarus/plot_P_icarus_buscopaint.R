
setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/2_sex_chr')

source('plot_buscopainter_functions.R')

complete_locations <- '../1_genomes/Polyommatus_icarus_complete_location.tsv'
prefix <- 'P. icarus'
index <- 'Polyommatus_icarus.seqlen.bed'
minimum <- 5
ncol = 1

location_set <- prepare_data_with_index(complete_locations, index)
locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(locations_filt)

locations_filt_Z_only <- locations_filt[locations_filt$query_chr == 'OW569320.1',]
Picarus_merian_differences_only_plot <- paint_merians_differences_only(spp_df, subset_merians, num_col, '', num_contigs)
paint_merians_all(locations_filt_Z_only, num_col, prefix, num_contigs)

ggsave(Picarus_merian_differences_only_plot, filename='Picarus_buscopaint_differences_only.png')
ggsave(Picarus_merian_differences_only_plot, filename='Picarus_buscopaint_differences_only.pdf')
MZ_locs <- locations_filt[(locations_filt$query_chr =="OW569320.1") & (locations_filt$status=="self"), ]
M24_locs <- locations_filt[(locations_filt$query_chr =="OW569320.1") & (locations_filt$status=="M24"), ] 

max(M24_locs$position) / 1000000 # 13.81665
min(MZ_locs$position) / 1000000 # 14.09043
# check to see if M24 portion is completely distinct from rest of MZ
# i.e. we can use a chr position cutoff to get M24-derived genes
ggplot(data=MZ_locs, aes(x=position, y=1.3)) + geom_point() + 
  geom_point(data=M24_locs, aes(x=position, y=1.4)) + ylim(1,2)

# so for P. icarus, we can define the M24 portion as from 1 to 13.81665 Mb


