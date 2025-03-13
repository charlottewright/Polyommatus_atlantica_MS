
library(dplyr)
library(ggplot2)
library(patchwork)
source('../2_sex_chr/plot_buscopainter_functions.R')
setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/5_telomeres/')

complete_locations <- 'Chrysodeixis_includens_complete_location.tsv'
prefix <- 'Chrysodeixis includens'

minimum <- 5
ncol = 1

location_set <- prepare_data(complete_locations)
locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(locations_filt)

locations_filt$assigned_chr_copy <- locations_filt$assigned_chr

Cincludens_merian_plot <- paint_merians_all(locations_filt, num_col, prefix, num_contigs)
Cincludens_merian_plot <- picarus_merian_plot + ggtitle('A') + theme(plot.title = element_text(face = "bold", hjust=0)) 

Cincludens_merian_plot

complete_locations <- 'Polyommatus_iphigenia_complete_location.tsv'
prefix <- 'Polyommatus iphigenia'

minimum <- 5
ncol = 1

location_set <- prepare_data(complete_locations)
locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(locations_filt)

locations_filt$assigned_chr_copy <- locations_filt$assigned_chr
Piphigenia_plot <- paint_merians_all(locations_filt, num_col, prefix, num_contigs)
Piphigenia_plot <- Piphigenia_plot + ggtitle('B') + theme(plot.title = element_text(face = "bold", hjust=0)) 



supfig_plot <- Cincludens_merian_plot + Piphigenia_plot + plot_layout(guides = 'collect')
ggsave(supfig_plot, filename='SupFig_C_includens_and_P_iphigenia_buscopaints.pdf',height=10, width=10)


