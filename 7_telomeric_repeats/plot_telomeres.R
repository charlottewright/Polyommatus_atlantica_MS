
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ape)
library(phytools)

# Aim: compare telomere prevalence and length in P. atlantica vs all other leps
setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/5_telomeres/')
chr_IDs <- read.csv('list_of_all_chr_IDs.tsv', sep=' ') # this gives the scaffold IDs for the chr (i.e. allows you to filter out non-chr shrapnel)
colnames(chr_IDs) <- c('chr', 'species')
head(chr_IDs)
P_atlantica_chr_lengths <- read.csv('Polyommatus_atlantica.chr_lengths.tsv',sep='\t', header=FALSE, col.names = c('chr','length'))

telos <- read.csv('telomeric_motifs_detected_by_TRF_533leps.080924.tsv', sep='\t', header=FALSE)[,c(1,2,3,4,16,18,19)]
colnames(telos) <- c('species', 'chr', 'start', 'stop', 'repeat', 'type', 'position')
telos$span <- telos$stop - telos$start
telos_filt <- telos %>% filter(chr %in% chr_IDs$chr)
head(telos_filt)
nrow(telos_filt) # 277,780 locations of telomeric sequence in chr
telos_end <- telos_filt %>% filter(type == 'end')
telos_middle <- telos_filt %>% filter(type == 'middle')
nrow(telos_end) # 214,388 locations of telomeric sequence are at ends of chr
nrow(telos_middle) # 63,392 locations of telomeric sequence are in middles of chr
mean(telos_middle$position) #mean position of middle telomere is 0.5!
length(unique(telos_middle$species)) # 530 species have >= 1 internal telomere lol

telos_middle_300bp_cutoff <- telos_middle %>% filter(span >= 300)
nrow(telos_middle_300bp_cutoff) # 2767
length(unique(telos_middle_300bp_cutoff$species)) # 134 species

# sum the number of ALL middle telomeres per species
species_count <- table(telos_middle$species)
species_count_df <- as.data.frame(species_count)
colnames(species_count_df) <- c("Species", "Count")
# make sure 'Count' is numeric
species_count_df$Count <- as.numeric(species_count_df$Count)

# sum the number of large (>300 bp) telomeres per species
species_count_300bp <- table(telos_middle_300bp_cutoff$species)
species_count_300bp_df <- as.data.frame(species_count_300bp)
colnames(species_count_300bp_df) <- c("Species", "Count")
# make sure 'Count' is numeric
species_count_300bp_df$Count <- as.numeric(species_count_300bp_df$Count)

# how many species have under 10 large ITS regions?
nrow(species_count_300bp_df %>% filter(Count <10))

# whats the proportion of internal telomeric arrays that are short?
short_telos_middle <- telos_middle %>% filter(span < 300)
(nrow(short_telos_middle) / nrow(telos_middle)) *100 # shows 96%!

# Plot count of ALL middle telomeres per species
all_telomere_count_plot <- ggplot(species_count_df, aes(x=Count)) + 
  ggtitle('A') + theme(plot.title = element_text(face = "bold")) +
  geom_histogram(bins=30,color = "black", fill = "black") + theme_bw() + labs(x="Number of internal telomeric arrays", y="Species count")

# Plot count of large middle telomeres per species
large_telomere_count_plot <- ggplot(species_count_300bp_df, aes(x=Count)) +
  ggtitle('B') + theme(plot.title = element_text(face = "bold"))+
  geom_histogram(bins=30,color = "black", fill = "black") + theme_bw() + labs(x="Number of large internal telomeric arrays", y="Species count")

sup_fig_count_dist <- all_telomere_count_plot + large_telomere_count_plot + plot_layout(nrow=2)
sup_fig_count_dist

ggsave(sup_fig_count_dist, filename='SupFig_histogram_of_number_telos_per_species.pdf')
ggsave(sup_fig_count_dist, filename='SupFig_histogram_of_number_telos_per_species.png')

# Plot a histogram of positions of ALL telomeres
hist_of_all_telos_pos <- ggplot(telos_filt, aes(x = position)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "black") +
  geom_vline(xintercept = 0.025, color = "grey", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.975, color = "grey", linetype = "dashed", size = 1) +
  labs(title = "Histogram of positions of telomeres", x = "Proportion of chr length", y = "Frequency") +
  theme_bw()

# Plot a histogram of the span of all internal telomeres
hist_of_span_of_internal_telos <- ggplot(telos_middle, aes(x = span)) +
  geom_histogram(color = "#7FC29B", fill = "#7FC29B") + # did have binwidth=100 when didn't use logscale
  labs(title = "Histogram of internal telomeric sequence span across all species", x = "Size of internal telomeric array (bases)", y = "Frequency") +
  theme_bw() + scale_x_log10()
# see that a few large ITS make x-axis huge but most ITS < 1000

p_atlantica_telos <- telos_filt %>% filter(species=="Polyommatus_atlantica")
p_atlantica_telos <- merge(p_atlantica_telos, P_atlantica_chr_lengths,by='chr')
p_atlantica_telos$prop_start <- p_atlantica_telos$start / p_atlantica_telos$length
p_atlantica_telos$prop_stop <- p_atlantica_telos$stop / p_atlantica_telos$length
p_atlantica_telos$prop_length <- p_atlantica_telos$stop / p_atlantica_telos$length

p_atlantica_middle_telos <- p_atlantica_telos %>% filter(type=="middle")
p_atlantica_end_telos <- p_atlantica_telos %>% filter(type=="end")
p_atlantica_end_1_telos <- p_atlantica_end_telos %>% filter(position < 0.025)
p_atlantica_end_2_telos <- p_atlantica_end_telos %>% filter(position > 0.975)

print('Number of chr ends with at least one telomeric array:')
number_ends_telos_1 <- length(unique(p_atlantica_end_1_telos$chr))
number_ends_telos_2 <- length(unique(p_atlantica_end_2_telos$chr))
print(number_ends_telos_1+number_ends_telos_2) # 247
print('Percentage of chr ends with an array:')
print(((number_ends_telos_1+number_ends_telos_2)/462)*100) # 94%
# total ends = 229+2 (for Ws) * 2 = 462

print('Number of internal telomeric arrays:')
nrow(p_atlantica_middle_telos) # 727 internal telomeric arrays

print('Number of chr with at least one ITS:')
print(length(unique(p_atlantica_middle_telos$chr)))

# lets count the number of internal telomeric arrays per chr
chr_count <- table(p_atlantica_middle_telos$chr)
chr_count_df <- as.data.frame(chr_count)
colnames(chr_count_df) <- c('chr', 'count')
max(chr_count_df$count) # 22
min(chr_count_df$count) # 1
nrow(chr_count_df %>% filter(count == 1))

test <- p_atlantica_telos %>% filter(chr %in% c('SUPER_1','SUPER_2','SUPER_3','SUPER_4', 'SUPER_5', 'SUPER_Z1', 'SUPER_Z2'))

theme_info <- theme(legend.position="right",
                    strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
                    panel.background = element_rect(fill = "white", colour = "white"), 
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line.x = element_line(color="black", size = 0.5),
                    # axis.text.x = element_text(size=15),
                    #  axis.title.x = element_text(size=15),
                    strip.text.y = element_text(angle=0),
                    strip.background = element_blank(),
                   # plot.title = element_text(hjust = 0.5, face="italic"),
                    plot.subtitle = element_text(hjust = 0.5))
# plot actual lengths of chr
ggplot(data=test) + 
  geom_rect(aes(xmin=0, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
  geom_rect(aes(xmin=start, xmax=stop, ymax=0, ymin =12), fill="blue") + 
  facet_wrap(chr ~., ncol=1, strip.position="right") + guides(scale="none") + 
  xlab("Position (Mb)") +
  scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
  scale_y_continuous(breaks=NULL) + theme_info

# plot as function of proportional length instead to allow autosomes+sex chr to be plotted together

plot_colours <- c('#BC4749','#7FC29B')
telomere_paint <- ggplot(data=test) + 
  geom_rect(aes(xmin=0, xmax=1, ymax=0, ymin =12), colour="black", fill="white") + 
  geom_rect(aes(xmin=prop_start-0.001, xmax=prop_stop+0.001, ymax=0, ymin =12, fill=type)) + 
  facet_wrap(chr ~., ncol=1, strip.position="right") + guides(scale="none") + 
  scale_fill_manual(values=plot_colours) +
  xlab("Position (Mb)") +  scale_y_continuous(breaks=NULL) + theme_info
telomere_paint

ggsave(telomere_paint, filename='telomere_locations_chr1_5_and_Zs.pdf')

summary_his_pos_and_size <- hist_of_all_telos_pos + theme(title=element_blank()) + hist_of_span_of_internal_telos + theme(title=element_blank()) + plot_layout(nrow=2)
summary_his_pos_and_size

fig_on_P_atlantica_plus_all_leps_summary <- telomere_paint + summary_his_pos_and_size + plot_annotation(title='A')
ggsave(fig_on_P_atlantica_plus_all_leps_summary, filename='telomere_paint_plus_summary_all_leps.pdf',width=10, height=6)
supfig_hist_of_all_telos_pos <- hist_of_all_telos_pos + theme(plot.title=element_blank())

ggsave(supfig_hist_of_all_telos_pos, filename='hist_of_all_telos_pos_across_leps.png', height=4, width=5)
write.table(p_atlantica_telos, quote = FALSE, row.names = FALSE, sep="\t",file = 'SupTable_locations_of_Patlantica_telomeric_arrays.tsv')


# lets use a cut-off of 300 bp to define short ITS
ggplot(short_telos_middle, aes(x = span)) +
  geom_histogram(binwidth = 10, color = "black", fill = "black") +
  labs(title = "Histogram of ITS span across all species", x = "Size of ITS (bases)", y = "Frequency") +
  theme_bw()

# Zoom in on the range of positions classified as "internal telomeres"
ggplot(telos_middle, aes(x = position)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "black") +
  labs(title = "Histogram of positions of internal telomeres", x = "Proportion of chr length", y = "Frequency") +
  theme_minimal()

mean_span_per_species <-telos_middle %>% group_by(species) %>% summarise(avg_value = mean(span))
mean_span_per_species <- data.frame(mean_span_per_species)
median_span_per_species <-telos_middle %>% group_by(species) %>% summarise(avg_value = median(span))
median_span_per_species <- data.frame(median_span_per_species)
mean_span_per_species <- mean_span_per_species[order(mean_span_per_species$avg_value, decreasing=TRUE),]
median_span_per_species <- median_span_per_species[order(median_span_per_species$avg_value, decreasing=TRUE),]
head(mean_span_per_species)
head(median_span_per_species)
# P. atlantica has the 6th largest mean span: 1032 bp
# P. atlantica has the 3rd largest median span 569 bp

# lets compare to sizes of telos at ends of chr
mean_end_span_per_species <-telos_end %>% group_by(species) %>% summarise(avg_value = mean(span))
mean_end_span_per_species <- data.frame(mean_end_span_per_species)
median_end_span_per_species <-telos_end %>% group_by(species) %>% summarise(avg_value = median(span))
median_end_span_per_species <- data.frame(median_end_span_per_species)
mean_end_span_per_species <- mean_end_span_per_species[order(mean_end_span_per_species$avg_value, decreasing=TRUE),]
median_end_span_per_species <- median_end_span_per_species[order(median_end_span_per_species$avg_value, decreasing=TRUE),]
head(mean_end_span_per_species)
head(mean_end_span_per_species)

total_span_per_species <- telos_middle %>%
  group_by(species) %>%
  summarize(sum_span = sum(span, na.rm = TRUE))
total_span_per_species <- data.frame(total_span_per_species)
# P. atlantica has a mean end span of 970 bp
# P. atlantica has a median end span of 534 bp
# lets combine stats to make a summary table as a suptable

summary_per_species <- merge(species_count_df, median_span_per_species, by.x='Species', by.y='species')
summary_per_species <- merge(summary_per_species, total_span_per_species, by.x='Species', by.y='species')

write.table(summary_per_species, quote = FALSE, row.names = FALSE, sep="\t",file = 'SupTable_summary_stats_telomeres_all_leps.tsv')


lycaenid_spp <- c('Polyommatus_atlantica', 'Polyommatus_icarus',
                  'Polyommatus_iphigenia','Lysandra_bellargus',
                  'Lysandra_coridon','Aricia_agestis','Aricia_artaxerxes',
                  'Cyaniris_semiargus','Plebejus_argus')

telos_middle_lycaenids <- telos_middle %>% filter(species %in% lycaenid_spp)
telos_middle_lycaenids$species <- as.factor(telos_middle_lycaenids$species)
telos_middle_lycaenids$species <- factor(telos_middle_lycaenids$species, levels = lycaenid_spp)

lycaenids_span_plot <- ggplot(telos_middle_lycaenids, aes(x=span)) + geom_histogram(binwidth =50, fill='#7FC29B') + 
  facet_wrap(~species, ncol=1) + theme_bw() + theme(strip.text = element_blank()) +
  xlim(0,4000) 

lycaenids_span_plot <- ggplot(telos_middle_lycaenids, aes(x=span)) + geom_histogram(fill='#7FC29B') + 
  facet_wrap(~species, ncol=1) + theme_bw() + theme(strip.text = element_blank()) +scale_x_log10()


lycaenids_span_plot
lycaenid_spp_with_long_ITS <-  c('Polyommatus_atlantica', 'Polyommatus_icarus',
                                 'Polyommatus_iphigenia','Lysandra_bellargus',
                                 'Lysandra_coridon')

telos_middle_lycaenids_long_ITS <- telos_middle_300bp_cutoff %>% filter(species %in% lycaenid_spp_with_long_ITS)
telos_middle_lycaenids_long_ITS$species <- as.factor(telos_middle_lycaenids_long_ITS$species)
telos_middle_lycaenids_long_ITS$species <- factor(telos_middle_lycaenids_long_ITS$species, levels = lycaenid_spp)

lycaenids_300bp_pos_plot <- ggplot(telos_middle_lycaenids_long_ITS, aes(x=position)) + 
  geom_histogram(fill='#7FC29B') + 
  facet_wrap(~species, ncol=1) + theme_bw() + theme(strip.text = element_blank()) 

combined_plot <- lycaenids_span_plot + lycaenids_300bp_pos_plot

ggsave(lycaenids_300bp_pos_plot, filename="internal_large_telomeres_pos_across_lycaenids.png")
ggsave(lycaenids_300bp_pos_plot, filename="internal_large_telomeres_pos_across_lycaenids.pdf")

ggsave(lycaenids_span_plot, filename="internal_telomeres_span_across_lycaenids.png")
ggsave(lycaenids_span_plot, filename="internal_telomeres_span_across_lycaenids.pdf")

ggsave(lycaenids_plot, filename="internal_telomeres_across_lycaenids2.png")
ggsave(lycaenids_plot, filename="internal_telomeres_across_lycaenids2.pdf")

ggsave(combined_plot, filename="internal_telomeres_across_lycaenids_and_pos_of_large_ITS.png")
ggsave(combined_plot, filename="internal_telomeres_across_lycaenids_and_pos_of_large_ITS.pdf", width=8, height=6)

tree <- read.tree('supermatrix.backtrans_nucleotides.MFP.rooted.txt')
tree <- drop.tip(tree, tip="Melitaea_cinxia")
tree <- drop.tip(tree, tip="Glaucopsyche_alexis")
tree <- drop.tip(tree, tip="Phengaris_arion")
tree <- drop.tip(tree, tip="Celastrina_argiolus")
tree <- drop.tip(tree, tip="Helleia_helle")
tree <- drop.tip(tree, tip="Lycaena_phlaeas")

cladogram_plot <- ggtree(tree, branch.length="none") + geom_tiplab() + ggplot2::xlim(0, 12) 
cladogram_plot

ggsave("lycaenids_subsetted_cladogram.pdf", plot = cladogram_plot, width = 5, height = 5)


# extra code to make figs
TRAS <- read.csv('P_atlantica.TRAS.bed', sep='\t', header=FALSE)
SART <- read.csv('P_atlantica.SART.bed', sep='\t', header=FALSE)
colnames(TRAS)<- c('chr', 'start', 'stop')
colnames(SART)<- c('chr', 'start', 'stop')

TRAS$type <- 'TRAS'
SART$type <- 'SART'
transposons <- rbind(TRAS, SART)
transposons <- merge(transposons, P_atlantica_chr_lengths,by='chr')

transposons$prop_start <- transposons$start / transposons$length
transposons$prop_stop <- transposons$stop / transposons$length

transposons_filt <- transposons %>% filter(chr %in% c('SUPER_1','SUPER_2','SUPER_3','SUPER_4', 'SUPER_5', 'SUPER_Z1', 'SUPER_Z2'))

# plot as function of proportional length instead to allow autosomes+sex chr to be plotted together
plot_colours <- c('#BC4749','#7FC29B', '#586994','#DB9D47')
telomere_plus_retrotransposons_paint <- ggplot(data=test) + 
  geom_rect(aes(xmin=0, xmax=1, ymax=0, ymin =6), colour="black", fill="white") + 
  geom_rect(aes(xmin=prop_start-0.001, xmax=prop_stop+0.001, ymax=0, ymin =6, fill=type)) + 
  facet_wrap(chr ~., ncol=1, strip.position="right") + guides(scale="none") + 
  scale_fill_manual(values=plot_colours) +
  xlab("Position (Mb)") +  scale_y_continuous(breaks=NULL) + theme_info +
  geom_rect(data=transposons_filt, aes(xmin=prop_start-0.001, xmax=prop_stop+0.001, ymax=6, ymin=9, fill=type))
telomere_plus_retrotransposons_paint

ggsave(telomere_plus_retrotransposons_paint, filename='telomere_locations_plus_SART_TRAS_chr1_5_and_Zs.pdf')

summary_num_per_spp_and_his_size <- all_telomere_count_plot + theme(title=element_blank()) + hist_of_span_of_internal_telos + theme(title=element_blank()) + plot_layout(nrow=2)
fig_on_P_atlantica_plus_all_leps_summary <- telomere_plus_retrotransposons_paint + summary_num_per_spp_and_his_size + plot_annotation(title='A')


# plot together
telos_middle_under_300bp_cutoff <- telos_middle %>% filter(span < 300)
species_count_small <- table(telos_middle_under_300bp_cutoff$species)
species_count_small <- as.data.frame(species_count_small)
colnames(species_count_small) <- c("Species", "Count")
species_count_small$type <- 'small'
species_count_300bp_df$type <- 'large'
test <- rbind(species_count_small,species_count_300bp_df)

col_palette <- c('black', 'darkgrey')
hist_num_ITS_per_spp <- ggplot(test, aes(x=Count, fill=type)) + geom_histogram() +  ggtitle('B') + theme(plot.title = element_text(face = "bold")) +
  theme_bw() + labs(x="Number of internal telomeric arrays", y="Species count") + scale_x_log10() +
  scale_fill_manual(values=col_palette)

B_and_C <- hist_num_ITS_per_spp +theme(plot.title=element_blank()) + hist_of_span_of_internal_telos + theme(plot.title=element_blank()) +plot_layout(nrow=2)
B_and_C
fig_on_P_atlantica_plus_all_leps_summary <- telomere_plus_retrotransposons_paint + B_and_C 
fig_on_P_atlantica_plus_all_leps_summary

ggsave(fig_on_P_atlantica_plus_all_leps_summary, filename='telomere_paint_with_retrotransposons_plus_summary_all_leps.log10scale.pdf',width=10, height=6)
