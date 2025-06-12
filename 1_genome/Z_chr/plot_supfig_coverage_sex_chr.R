
# Aim: Plot PacBio coverage per scaffold to identify sex-linked scaffolds
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)

busco_paint_one_chr <- function(spp_df,custom_chr_order = NULL){
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
                    'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18',
                    'M19', 'M20', 'M21', 'M22', 'M23', 'M24', 'M25', 'M26',
                    'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  
  # Ensure consistent factor levels for assigned_chr
  spp_df$assigned_chr <- factor(spp_df$assigned_chr, levels = merian_order)
  
  # Consistent color palette
  palette_colors <- setNames(append(hue_pal()(32), 'grey'), merian_order)
  
  # lets edit the palette
  palette_colors["MZ"] <- "#227c9d"
  palette_colors["M16"] <- "#ffcb77"
  palette_colors["M21"] <- "#17c3b2"
  palette_colors["M24"] <- "#fe6d73"
  palette_colors["M17"] <- "#399E5A"
  palette_colors["M20"] <- "#A7BBEC"
  # Consistent chromosome order
  # Determine chromosome order
  if (!is.null(custom_chr_order)) {
    chr_levels <- custom_chr_order
  } else {
    chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE)
    chr_levels <- chr_levels$query_chr
  }
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting
  print(chr_levels)
  the_plot <- ggplot(data = spp_df) + 
    scale_colour_manual(values = palette_colors, aesthetics = c("colour", "fill")) +
    geom_rect(aes(xmin = start, xmax = length, ymax = 0, ymin = 12), colour = "black", fill = "white") + 
    geom_rect(aes(xmin = position - 2e4, xmax = position + 2e4, ymax = 0, ymin = 12, fill = assigned_chr)) +
    facet_wrap(query_chr_f ~., nrow=1,scales="free") + guides(scale="none") + 
    theme(legend.position = "none",
          strip.text = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "white"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.text.x = element_text(size = 15,color='black'),
          axis.text.y = element_text(color='white'),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 15),
          strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    xlab('Position (Mb)') +
    scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.005, 1)) +
    scale_y_continuous(breaks = NULL)
  
  return(the_plot)
}



setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/1_genomes')
coverage_df <- read.csv('../2_sex_chr/Polyommatus_atlantica.100kb.mosdepth.regions.bed', sep='\t', header=FALSE)
repeats_df <- read.csv('Polyommatus_atlantica.filteredRepeats.per_scaffold.tsv', sep='\t', header=FALSE)[,c(1,3,7)]
colnames(coverage_df) <- c('scaff', 'start', 'end', 'coverage')
colnames(repeats_df) <- c('scaff', 'length', 'proportion')

#---------------------------------------------------------------------------------
# Fig A&B: First lets plot average coverage and repeat density per scaffold in P. atlantica
#---------------------------------------------------------------------------------
# calculate average coverage per scaffold
average_coverage <- coverage_df %>%
  group_by(scaff) %>%
  summarise(avg_coverage = mean(coverage))

# get size of each scaffold
scaff_size <- coverage_df %>%
  group_by(scaff) %>%
  summarise(length = max(end))

average_coverage_df <- merge(scaff_size, average_coverage, by='scaff')
average_coverage_df$length_mb <- average_coverage_df$length / 1000000
repeats_df$length_mb <- repeats_df$length / 1000000

# smallest chr is SUPER_227 which is 1226658 (1.2 Mb) - so lets filter for >1 Mb to remove shrapnel

average_coverage_df_filt <- average_coverage_df[average_coverage_df$length_mb >= 1.2,]
repeats_df_filt <- repeats_df[repeats_df$length_mb >= 1.2,]

#average_coverage_sex_chr <- average_coverage_df_filt[average_coverage_df_filt$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]
average_coverage_Z_chr <- average_coverage_df_filt[average_coverage_df_filt$scaff %in% c('SUPER_Z1', 'SUPER_Z2'),]
average_coverage_W_chr <- average_coverage_df_filt[average_coverage_df_filt$scaff %in% c('SUPER_W1', 'SUPER_W2'),]

average_coverage_autosome <- average_coverage_df_filt[!average_coverage_df_filt$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]
#average_coverage_sex_chr$type <- 'sex_chr'
average_coverage_autosome$Type <- 'Autosome'
average_coverage_Z_chr$Type <- 'Z'
average_coverage_W_chr$Type <- 'W'
average_coverage_filt_type_df <- rbind(average_coverage_autosome, average_coverage_Z_chr, average_coverage_W_chr)

# Do same for repeat density
repeats_Z_chr <- repeats_df_filt[repeats_df_filt$scaff %in% c('SUPER_Z1', 'SUPER_Z2'),]
repeats_W_chr <- repeats_df_filt[repeats_df_filt$scaff %in% c('SUPER_W1', 'SUPER_W2'),]
repeats_autosomes <- repeats_df_filt[!repeats_df_filt$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]
repeats_autosomes$Type <- 'Autosome'
repeats_Z_chr$Type <- 'Z'
repeats_W_chr$Type <- 'W'

repeats_filt_type_df <- rbind(repeats_autosomes, repeats_Z_chr, repeats_W_chr)
repeats_filt_type_df$repeat_per <- repeats_filt_type_df$proportion * 100

# Make plots for P. atlantica
cov_plot <- ggplot(average_coverage_filt_type_df, aes(x=length_mb, y=avg_coverage, color=Type)) + geom_point(size=2) + 
  scale_color_manual(values = c("W" = "red", "Z" = "blue", "Autosome" = "grey")) +
  labs(x="Length (Mb)", y = "Average coverage") + 
  theme_bw() + ggtitle('A') + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

repeats_plot <- ggplot(repeats_filt_type_df, aes(x=length_mb, y=repeat_per, color=Type)) + geom_point(size=2) + 
  scale_color_manual(values = c("W" = "red", "Z" = "blue", "Autosome" = "grey")) +
  labs(x="Length (Mb)", y = "Repeat density (%)") + 
  theme_bw() +  ggtitle('B') + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

cov_and_repeats_plot <- cov_plot + repeats_plot + plot_layout(guides = 'collect')
cov_and_repeats_plot

#---------------------------------------------------------------------------------
# FigC: Now lets plot pacbio coverage along chr of P. atlantica (use without Ws)
#---------------------------------------------------------------------------------
source('../functions_busco_painter.R')
#coverage_df <- read.csv('../2_sex_chr/Polyommatus_atlantica.100kb.mosdepth.regions.bed', sep='\t', header=FALSE)
coverage_noWs_df <- read.csv('../2_sex_chr/Polyommatus_atlantica.without_Ws.100kb.mosdepth.regions.bed', sep='\t', header=FALSE)
colnames(coverage_noWs_df) <- c('scaff', 'start', 'end', 'coverage')

# calculate average coverage per scaffold
average_coverage_no_Ws <- coverage_noWs_df %>%
  group_by(scaff) %>%
  summarise(avg_coverage = mean(coverage))

average_coverage_no_Ws_df <- merge(scaff_size, average_coverage_no_Ws, by='scaff')
autosomes_coverage_no_Ws_df <- coverage_noWs_df[!coverage_noWs_df$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]
coverage_noWs_df$Female <- coverage_noWs_df$coverage / median(autosomes_coverage_no_Ws_df$coverage)

subsetted_coverage_no_Ws_df <- coverage_noWs_df[coverage_noWs_df$scaff %in% c('SUPER_Z1', 'SUPER_Z2','SUPER_1'),]
subsetted_coverage_no_Ws_df$midpos <- (subsetted_coverage_no_Ws_df$start + subsetted_coverage_no_Ws_df$end)/2

sex_chr_windows_P_atlantica_noWs_plot <- ggplot(subsetted_coverage_no_Ws_df, aes(x=midpos/1000000, y=Female)) + geom_line(color='#337CA0') + 
  facet_wrap(~scaff, nrow=1, scales='free_x') + ylim(0,2.5) + labs(x='Position (Mb)', y=expression("Normalised coverage (" * italic("N") * ")")) +  
  theme_bw() + ggtitle(bquote(bold('C') ~ italic('    P. atlantica'))) + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15, color='black')) +  theme(axis.title = element_text(size=15, color='black')) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


complete_locations <- '../1_genomes/Polyommatus_atlantica_complete_location.tsv'
duplicate_locations <- '../1_genomes/Polyommatus_atlantica_duplicated_location.tsv'

prefix <- 'P. atlantica'

index <- '../1_genomes/Polyommatus_atlantica.fa.fai'
minimum <- 5
ncol = 1

location_set <- prepare_data_with_index(complete_locations, index)
dups_location_set <- prepare_data_with_index(duplicate_locations, index)
location_set <- rbind(location_set, dups_location_set)
Patlantica_locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(Patlantica_locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(Patlantica_locations_filt)

Patlantica_locations_filt$assigned_chr_copy <- Patlantica_locations_filt$assigned_chr


# pick the atlantica chr to make busco plots of
Patlantica_locations_filt$position_flipped <- Patlantica_locations_filt$length - Patlantica_locations_filt$position
Patlantica_target_chr <- Patlantica_locations_filt[Patlantica_locations_filt$query_chr %in% c('SUPER_1','SUPER_Z1','SUPER_Z2'),]

custom_chr_order <- c('SUPER_1','SUPER_Z1','SUPER_Z2')
Patlantica_target_chr_plot <- busco_paint_one_chr(Patlantica_target_chr,custom_chr_order)
Patlantica_target_chr_plot <- Patlantica_target_chr_plot + theme(plot.title = element_text(face = "bold",size=15))


panel_C <- sex_chr_windows_P_atlantica_noWs_plot + Patlantica_target_chr_plot + 
  plot_layout(nrow=2, heights = c(2,1)) + theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

#-----------------------------------
# Fig D: Lets plot P. icarus
#-----------------------------------
# Make plots for P. icarus
male_windows <- read.csv('../2_sex_chr/Polyommatus_icarus.male.mapped_to_male.100kb.mosdepth.regions.bed',sep='\t',header=FALSE)
female_windows <- read.csv('../2_sex_chr/Polyommatus_icarus.female.mapped_to_male.100kb.mosdepth.regions.bed',sep='\t',header=FALSE)
colnames(male_windows) <- c('chr', 'start','end','male_cov')
colnames(female_windows) <- c('chr', 'start','end','female_cov')
male_windows$midpos <- (male_windows$start + male_windows$end) / 2
female_windows$midpos <- (female_windows$start + female_windows$end) / 2

# calculate median coverage for autosomes
male_autosomes <- male_windows[male_windows$chr != 'OW569320.1',]
female_autosomes <- female_windows[female_windows$chr != 'OW569320.1',]

# now normalise raw coverage values per median coverage per chr
male_windows$Male <- male_windows$male_cov / median(male_autosomes$male_cov)
female_windows$Female <- female_windows$female_cov / median(female_autosomes$female_cov)

windows_merged <- merge(male_windows, female_windows, by=c('chr', 'start','end','midpos'))

windows_merged_long <- windows_merged %>%
  pivot_longer(cols = c(Male, Female), 
               names_to = "sex", 
               values_to = "cov_norm")

# Z chr in P. icarus is OW569320.1 (M24+MZ)
# lets check OW569334.1 (M16) which is Z2 in P. atlantica
# show an autosome for comparison: M20 = OW569328.1
subsetted_windows_merged_long <- windows_merged_long[windows_merged_long$chr %in% c('OW569320.1', 'OW569334.1','OW569328.1'),]
chr_order <- c('OW569328.1', 'OW569320.1','OW569334.1')
subsetted_windows_merged_long$chr_f <- factor(subsetted_windows_merged_long$chr, levels=chr_order)

P_icarus_chr_lengths <- subsetted_windows_merged_long %>% group_by(chr) %>% summarise(max(end))
colnames(P_icarus_chr_lengths) <- c('chr', 'chr_length')

subsetted_windows_merged_long <- merge(subsetted_windows_merged_long, P_icarus_chr_lengths, by='chr')

subsetted_windows_merged_long$midpos_flipped <- subsetted_windows_merged_long$chr_length - subsetted_windows_merged_long$midpos

sex_chr_Picarus_flipped_plot <- ggplot(subsetted_windows_merged_long, aes(x=midpos_flipped/1000000, y=cov_norm, color=sex)) + geom_line() + 
  facet_wrap(~chr_f, nrow=1, scales='free_x') + ylim(0,2.5) + labs(x='Position (Mb)', y=expression("Normalised coverage (" * italic("N") * ")")) +  
  scale_color_manual(values=c('#337CA0', '#0D090A')) +
  theme_bw() + ggtitle(bquote(bold('D') ~ italic('    P. icarus'))) + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15, color='black')) +  theme(axis.title = element_text(size=15,color='black')) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


# pick the P.icarus chr to make busco plots of
source('../functions_busco_painter.R')

complete_locations <- '../2_sex_chr/Polyommatus_icarus_complete_location.tsv'
prefix <- 'P. icarus'

index <- '../2_sex_chr/Polyommatus_icarus.fa.fai'
minimum <- 5
ncol = 1

location_set <- prepare_data_with_index(complete_locations, index)
locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(locations_filt)

locations_filt$assigned_chr_copy <- locations_filt$assigned_chr
P_icarus_chr_lengths <- windows_merged_long %>% group_by(chr) %>% summarise(max(end))
locations_filt <- merge(locations_filt, P_icarus_chr_lengths, by.x='query_chr',by.y='chr')
locations_filt$position_flipped <- locations_filt$length - locations_filt$position
P_icarus_OW569320_target_chr <- locations_filt[locations_filt$query_chr %in% c('OW569320.1','OW569328.1','OW569334.1'),]
P_icarus_OW569320_target_chr <- P_icarus_OW569320_target_chr %>%
  mutate(position = if_else(query_chr == "OW569320.1", position_flipped, position))


custom_chr_order <- c('OW569328.1','OW569320.1','OW569334.1')
P_icarus_OW569320_target_chr_plot <- busco_paint_one_chr(P_icarus_OW569320_target_chr,custom_chr_order)
P_icarus_OW569320_target_chr_plot <- P_icarus_OW569320_target_chr_plot + theme(plot.title = element_text(face = "bold",size=15))


panel_D <- sex_chr_Picarus_flipped_plot + P_icarus_OW569320_target_chr_plot + plot_layout(nrow=2, heights = c(2,1)) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))


#----------------------------------------------------------
# Fig E: Busco paint for all chromosomes in P. icarus
#----------------------------------------------------------

paint_merians_all <- function(spp_df, num_col, title, karyotype){
  #colour_palette <- append(hue_pal()(32), 'grey')
  colour_palette <- c('#227c9d', '#ffcb77',  '#17c3b2', '#fe6d73','grey')
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self','Other')
  spp_df$assigned_chr_f =factor(spp_df$assigned_chr, levels=merian_order)
  chr_levels = unique(spp_df$query_chr)
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting
  #sub_title <- paste("n contigs =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=colour_palette, aesthetics=c("colour", "fill")) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=assigned_chr_f)) +
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + 
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title)  +
    guides(fill=guide_legend("Merian element"), color = "none") +
    busco_paint_theme
  return(the_plot)
}


locations_filt_simplified <- locations_filt
locations_filt_simplified$assigned_chr[locations_filt_simplified$assigned_chr != 'M24' & locations_filt_simplified$assigned_chr != 'MZ'  & locations_filt_simplified$assigned_chr != 'M21'  & locations_filt_simplified$assigned_chr != 'M16'] <- "Other"

locations_filt_simplified <- locations_filt_simplified %>%
  mutate(position = if_else(query_chr == "OW569320.1", position_flipped, position))

prefix <- ''
picarus_merian_plot <- paint_merians_all(locations_filt_simplified, num_col, prefix, num_contigs) + 
  ggtitle('E') + theme(plot.title = element_text(face = "bold", hjust=0,size=15)) +
  theme(axis.text = element_text(color='black'))

#----------------------------------------------------------
# Fig F: Oxford plot for Z chr of P. atlantica vs P. icarus
#----------------------------------------------------------

source('../functions_oxford_plots.R')
Merian_assignments_ref <- read_busco('../2_sex_chr/Merian_elements_full_table.tsv')
Merian_assignments_ref <- select(Merian_assignments_ref, -c(start, end))
colnames(Merian_assignments_ref) <- c('busco_id', 'Merian')

species_1 <- read_duplicated_and_complete_buscos('Polyommatus_atlantica.tsv') # species on Y-axis
species_2 <- read_duplicated_and_complete_buscos('Polyommatus_icarus.tsv') # species on X-axis

species_1 <- species_1[species_1$Sequence %in% c('SUPER_Z1'),] # species on Y-axis
species_2 <- species_2[species_2$Sequence %in% c('OW569320.1'),] # species on X-axis

# lets flip P. icarus to make comparable to P. atlantica
species_2 <- species_2 %>% group_by(Sequence) %>% 
  mutate(seq_length = max(start, na.rm = TRUE)) %>%
  ungroup()

species_2$start <- species_2$seq_length - species_2$start
species_2$end <- species_2$seq_length - species_2$end

species_1_name <- 'Chromosome length (Mb)'
species_2_name <- 'Chromosome length (Mb)'

species_1 <- scale_lengths(species_1)
species_2 <- scale_lengths(species_2)

species_1_vs_species_2 <- make_pairwise_dataframe(species_1, species_2)
species_1_vs_species_2 <- filter_buscos_x(species_1_vs_species_2,2) %>% filter_buscos_y(1)
species_1_vs_species_2 <- merge(species_1_vs_species_2, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
species_1_vs_species_2$Merian <- factor(species_1_vs_species_2$Merian, levels = merian_order)

# try filtering by minimum merian element
species_1_vs_species_2 <- species_1_vs_species_2 %>% 
  group_by(Merian, Sequence.x) %>%
  mutate(nGenes = n()) %>%
  ungroup() %>%
  filter(nGenes >1)

# Make columns to define limits of each facet
species_1_vs_species_2 <- species_1_vs_species_2 %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
species_1_vs_species_2 <- species_1_vs_species_2 %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
species_1_vs_species_2$sequence_y_min <- 0
species_1_vs_species_2$sequence_x_min <- 0

# Default order of chr is by size 
sequence_x_order <- species_1_vs_species_2 %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- species_1_vs_species_2 %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()
sequence_x_order <- rev(sequence_x_order)

species_1_vs_species_2$Sequence.x_f <- factor(species_1_vs_species_2$Sequence.x, levels = sequence_x_order)
species_1_vs_species_2$Sequence.y_f <- factor(species_1_vs_species_2$Sequence.y, levels = sequence_y_order)

col_palette <- c('#227c9d', '#ffcb77', '#17c3b2','#fe6d73')
# correct function
make_oxford_plot <- function(df){
  the_plot <- ggplot() +  # scales=free, space=free is what you want to get boxes scaled by length :) 
    facet_grid(Sequence.x_f~Sequence.y_f, space = "free", scales="free")  + 
    geom_point(data = df, aes(x = midpos_y, y = midpos_x, color=Merian), size=0.1) + # was Sequence.x_f 
    theme_bw() + 
    theme(legend.position = "none", axis.ticks = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_rect(color = "#999999", size=0.3)) + 
    theme(panel.spacing = unit(0, "lines")) +
    theme(strip.background = element_blank()) +
    theme(strip.text.y = element_text(angle = 0,size=5)) +
    theme(strip.text.x = element_text(angle = 90,size=5)) +
    geom_point(data=df, aes(x=sequence_y_min, y=sequence_x_min), color="white") + # needed to scale each box
    geom_point(data=df, aes(x=sequence_y_max, y=sequence_x_max), color="white")     # needed to scale each box
  return(the_plot)
}

P_icarus_Z_vs_P_atlantica_Z <- make_oxford_plot(species_1_vs_species_2) + 
  xlab(species_2_name) +   ylab(species_1_name) +
  theme(legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=2), 
                             title="Merian element", nrow=2,byrow=TRUE)) +
  scale_color_manual(values=col_palette) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

P_icarus_Z_vs_P_atlantica_Z





#----------------------------------------------------------
# Final figure
#----------------------------------------------------------
layout <- "
AB
AB
AB
CC
CC
CC
DD
EE
EE
EE
FF
"

supfig_pt1 <- cov_plot + repeats_plot + sex_chr_windows_P_atlantica_noWs_plot + Patlantica_target_chr_plot + 
  sex_chr_Picarus_flipped_plot + P_icarus_OW569320_target_chr_plot + 
  plot_layout(guides = 'collect', design=layout)

supfig_pt2 <- picarus_merian_plot + theme(legend.position = "none") + P_icarus_Z_vs_P_atlantica_Z + plot_layout(guides = 'collect', widths=c(2,1)) +
  theme(legend.position = "none")


ggsave(supfig_pt1, filename='sex_chr_coverage_sup_fig_pt1.300525.pdf')
ggsave(supfig_pt2, filename='sex_chr_coverage_sup_fig_pt2.300525.pdf')
