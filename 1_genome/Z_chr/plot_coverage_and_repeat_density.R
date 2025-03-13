
# Aim: Plot PacBio coverage per scaffold to identify sex-linked scaffolds

library(dplyr)
library(ggplot2)
library(patchwork)

setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/1_genomes')
coverage_df <- read.csv('../2_sex_chr/Polyommatus_atlantica.100kb.mosdepth.regions.bed', sep='\t', header=FALSE)
repeats_df <- read.csv('Polyommatus_atlantica.filteredRepeats.per_scaffold.tsv', sep='\t', header=FALSE)[,c(1,3,7)]
colnames(coverage_df) <- c('scaff', 'start', 'end', 'coverage')
colnames(repeats_df) <- c('scaff', 'length', 'proportion')

# calculate average coverage per scaffold
average_coverage <- coverage_df %>%
  group_by(scaff) %>%
  summarise(avg_coverage = mean(coverage))

# get size of each scaffold
scaff_size <- coverage_df %>%
  group_by(scaff) %>%
  summarise(length = max(end))

head(scaff_size)
head(average_coverage)

average_coverage_df <- merge(scaff_size, average_coverage, by='scaff')
average_coverage_df$length_mb <- average_coverage_df$length / 1000000
repeats_df$length_mb <- repeats_df$length / 1000000

# smallest chr is SUPER_227 which is 1226658 (1.2 Mb)
# so lets filter for >1 Mb to remove shrapnel

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

ggsave(cov_and_repeats_plot, filename='pacbio_coverage_and_repeat_density_per_chr.pdf', device = 'pdf', width = 8, height = 4)
ggsave(cov_and_repeats_plot, filename='pacbio_coverage_and_repeat_density_per_chr.png', device = 'png', width = 8, height = 4)

coverage_df$midpos <- (coverage_df$start + coverage_df$end) / 2

autosomes_coverage_df <- coverage_df[!coverage_df$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]
coverage_df$Female <- coverage_df$coverage / median(autosomes_coverage_df$coverage)

subsetted_coverage_df <- coverage_df[coverage_df$scaff %in% c('SUPER_Z1', 'SUPER_Z2','SUPER_1'),]

sex_chr_windows_P_atlantica_plot <- ggplot(subsetted_coverage_df, aes(x=midpos/1000000, y=Female)) + geom_line(color='#337CA0') + 
  facet_wrap(~scaff, nrow=1, scales='free') + ylim(0,2.5) + labs(x='Position (Mb)', y='Normalised coverage (N)') + 
  theme_bw() + ggtitle('C') + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

sex_chr_windows_P_atlantica_plot

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

# see that the first half (up to 14 Mb) has similar coverage in M & F, think as W-linked sequences are also mapping 
sex_chr_Picarus_plot <- ggplot(subsetted_windows_merged_long, aes(x=midpos/1000000, y=cov_norm, color=sex)) + geom_line() + 
  facet_wrap(~chr_f, nrow=1, scales='free') + ylim(0,2.5) + labs(x='Position (Mb)', y=expression("Normalised coverage (" * italic("N") * ")")) +  
  scale_color_manual(values=c('#337CA0', '#0D090A')) + 
  theme_bw() + ggtitle('D') + theme(plot.title = element_text(face = "bold",size=15)) +
  theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))

# now lets combine coverage with buscopaint of Merians that are sex-linked in P. atlantica

setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/2_sex_chr')

source('plot_buscopainter_functions.R')

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

complete_locations <- 'Polyommatus_icarus_complete_location.tsv'
duplicated_locations <- 'Polyommatus_icarus_duplicated_location.tsv'
#complete_locations <- 'Polyommatus_icarus.female_complete_location.tsv'
#duplicated_locations <- 'Polyommatus_icarus.female_duplicated_location.tsv'
prefix <- 'P. icarus'

index <- 'Polyommatus_icarus.fa.fai'
minimum <- 5
ncol = 1

location_set <- prepare_data_with_index(complete_locations, index)
#location_set_dups <- prepare_data_with_index(duplicated_locations, index)
#location_set <- rbind(location_set, location_set_dups)
locations_filt <- filter_buscos(location_set, minimum) 
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_col <- 1
subset_merians <- set_merian_colour_mapping(locations_filt)

locations_filt$assigned_chr_copy <- locations_filt$assigned_chr
locations_filt$assigned_chr[locations_filt$assigned_chr != 'M24' & locations_filt$assigned_chr != 'MZ'  & locations_filt$assigned_chr != 'M21'  & locations_filt$assigned_chr != 'M16'] <- "Other"

#plot all
prefix <- ''
picarus_merian_plot <- paint_merians_all(locations_filt, num_col, prefix, num_contigs)
picarus_merian_plot <- picarus_merian_plot + ggtitle('E') + theme(plot.title = element_text(face = "bold", hjust=0,size=15)) 
P_icarus_cov_plot <-  P_icarus_cov_plot + theme(axis.text = element_text(size=15)) +  theme(axis.title = element_text(size=15))


layout <- "
AB
CC
DD
EE
EE
"


supfig_plot <- cov_plot + repeats_plot + sex_chr_windows_P_atlantica_plot + sex_chr_Picarus_plot + picarus_merian_plot + plot_layout(guides = 'collect', design=layout)
supfig_plot



# plot full SupFig of just Merians associated with sex chr in P. atlantica (complete and duplicated loci) and coverage plot 
ggsave(supfig_plot, filename='P_atlantica_and_P_icarus_sex_chr_identifiation_fig_using_M_vs_F_pacbio.260225.pdf',height=18, width=12)


