library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/1_genomes/')

# functions
get_merian_assignments_for_chr_with_few_markers <- function(df, max_markers){
  species_chr2merian <- df[df$status == 'self',]
  species_chr2merian <- species_chr2merian[,c(2,4)]
  colnames(species_chr2merian) <- c('scaff', 'Merian_element')
  species_chr2merian_counts <- species_chr2merian %>%
    group_by(scaff, Merian_element) %>%
    summarise(count = n(), .groups = 'drop')
  species_chr2merian_counts <- species_chr2merian_counts[species_chr2merian_counts$count < max_markers,]
  return(species_chr2merian_counts)
}

# Plot P. atlantica & P. icarus busco paints

# Merian locations
Picarus <- read.csv('../2_sex_chr/Polyommatus_icarus_complete_location.tsv', sep='\t')
Patlantica <- read.csv('Polyommatus_atlantica_complete_location.tsv', sep='\t')
#Csemiargus_assignments <- read.csv('Cyaniris_semiargus_chromosome_assignments.tsv', sep='\t')
Picarus_chr_assignments <- read.csv('Polyommatus_icarus.window_17_chromosome_assignments.tsv', sep='\t')
Patlantica_chr_assignments <- read.csv('Polyommatus_atlantica.window_3_chromosome_assignments.tsv', sep='\t')

# Chr lengths
Cyaniris_lengths <- read.csv('Cyaniris_semiargus.chr_lengths.bed', sep='\t', header=FALSE)[,c(1,3)]
Picarus_lengths <- read.csv('Polyommatus_icarus.chr_lengths.bed', sep='\t', header=FALSE)[,c(1,3)]

# Repeat density
repeats_df <- read.csv('Polyommatus_atlantica.filteredRepeats.per_scaffold.tsv', sep='\t', header=FALSE)[,c(1,3,7)]
genes_df <- read.csv('Polyommatus_atlantica.braker1_and_braker2.TSEBRA.gene_density_per_scaff.tsv',sep='\t', header=FALSE)[,c(1,3,4,7)]

# Merian assignments from LFFF with a window_size of 3
colnames(Patlantica_chr_assignments) <- c('scaff', 'status', 'Merian_element')

# get Merian assignments that weren't detected with a cut-off with a window size of 3
under_three_markers_Patlantica_chr2merian <- get_merian_assignments_for_chr_with_few_markers(Patlantica, 3) # DOESNT CONSIDER FUSIONS!
under_three_markers_Patlantica_chr2merian <-  under_three_markers_Patlantica_chr2merian %>%
  filter(grepl("SUPER", scaff)) %>%
  filter(!grepl("SUPER_W2", scaff)) %>%
  filter(!grepl("unloc", scaff))
under_three_markers_Patlantica_chr2merian <- under_three_markers_Patlantica_chr2merian[,c(1,2)]
Patlantica_assignments_long <- separate_rows(Patlantica_chr_assignments, Merian_element, sep = ",")[,c(1,3)] # make two rows per fusion chr
Patlantica_assignments_long <- rbind(Patlantica_assignments_long, under_three_markers_Patlantica_chr2merian)
Patlantica_assignments_long <- unique(Patlantica_assignments_long)

colnames(Cyaniris_lengths) <- c('scaff', 'length')
colnames(Picarus_lengths) <- c('scaff', 'length')
colnames(repeats_df) <- c('scaff', 'length', 'proportion')
colnames(genes_df) <- c('scaff', 'length', 'gene_count', 'proportion')

#Csemiargus_assignments_long <- separate_rows(Csemiargus_assignments, assigned_ref_chr, sep = ",")[,c(1,3)] # make two rows per fusion chr
#colnames(Csemiargus_assignments_long) <- c('scaff', 'Merian_element')

Picarus_assignments_long <- separate_rows(Picarus_chr_assignments, assigned_ref_chr, sep = ",")[,c(1,3)] # make two rows per fusion chr
colnames(Picarus_assignments_long) <- c('scaff', 'Merian_element')

repeats_df$length_mb <- repeats_df$length / 1000000
repeats_df$repeat_per <- repeats_df$proportion * 100
genes_df$length_mb <- genes_df$length / 1000000

# filter out unlocalised scaffolds
repeats_df <- repeats_df[repeats_df$length_mb >= 1.2,]
genes_df <- genes_df[genes_df$length_mb >= 1.2,]

Picarus_merian2length <- merge(Picarus_assignments_long, Picarus_lengths, by="scaff")
colnames(Picarus_merian2length) <- c('scaff_Picarus', 'Merian_element', 'length_Picarus')

#Csemiargus_merian2length <- merge(Csemiargus_assignments_long, Cyaniris_lengths, by="scaff")
#colnames(Csemiargus_merian2length) <- c('scaff_Csemiargus', 'Merian_element', 'length_Csemiargus')

repeats_df <- merge(repeats_df, Patlantica_assignments_long, by='scaff')
#repeats_df <- merge(repeats_df, Csemiargus_merian2length, by='Merian_element')
temp <- merge(repeats_df, Picarus_merian2length, by='Merian_element')

repeats_df <- merge(repeats_df, Picarus_merian2length, by='Merian_element')

merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
repeats_df$Merian_element_f =factor(repeats_df$Merian_element, levels=merian_order)
repeats_df_autosomes <- repeats_df[!repeats_df$scaff %in% c('SUPER_W1', 'SUPER_W2', 'SUPER_Z1', 'SUPER_Z2'),]

# count number of chr per merian
merian_counts <- repeats_df_autosomes %>%
  group_by(Merian_element) %>%
  summarise(count = n())

cat('Average number of chr per Merian:', mean(merian_counts$count))
cat('Min number of chr per Merian:', min(merian_counts$count))
cat('Max number of chr per Merian:', max(merian_counts$count))


# check for correlation between repeat density in P. atlantica and length of the chr in P. atlantica
length_vs_repeat_corr_test <- cor.test(repeats_df_autosomes$length, repeats_df_autosomes$repeat_per, method = "spearman")
length_vs_repeat_corr <- length_vs_repeat_corr_test$estimate # R: -0.08170576 
length_vs_repeat_pval <- length_vs_repeat_corr_test$p.value # p-val: 0.2120551
cat('R:', length_vs_repeat_corr, 'P-val:', length_vs_repeat_pval, '\n')
# conclude: no correlation!

# check for correlation between repeat density in P. atlantica and length of the chr in C. semiargus
#repeats_df_autosomes$length_Csemiargus <- repeats_df_autosomes$length_Csemiargus / 1000000
#length_vs_repeat_corr_test <- cor.test(repeats_df_autosomes$length_Csemiargus, repeats_df_autosomes$repeat_per, method = "spearman")
#length_vs_repeat_corr <- length_vs_repeat_corr_test$estimate # -0.2746755
#length_vs_repeat_pval <- length_vs_repeat_corr_test$p.value #  1.948841e-05 
#cat('R:', length_vs_repeat_corr, 'P-val:', length_vs_repeat_pval, '\n')
# conclude: strong correlation!

# check for correlation between repeat density in P. atlantica and length of the chr in P. icarus
repeats_df_autosomes$length_Picarus <- repeats_df_autosomes$length_Picarus / 1000000
length_vs_repeat_corr_test <- cor.test(repeats_df_autosomes$length_Picarus, repeats_df_autosomes$repeat_per, method = "spearman")
length_vs_repeat_corr <- length_vs_repeat_corr_test$estimate # R: -0.2651694
length_vs_repeat_pval <- length_vs_repeat_corr_test$p.value #p-val 3.834381e-05
cat('R:', length_vs_repeat_corr, 'P-val:', length_vs_repeat_pval, '\n')
# conclude: strong correlation!


# plot repeat density vs chr length in C. semiargus
repeats_vs_length_in_Csemiargus <- 
  ggplot(repeats_df_autosomes, aes(x=length_Csemiargus, y=repeat_per)) + 
  geom_point(size = 2, aes(color=factor(Merian_element_f))) +
  labs(x="Chr length in C. semiargus (Mb)", y = "Repeat density (%)") + 
  theme_classic() +  ggtitle('D') + theme(plot.title = element_text(face = "bold")) +
  geom_smooth(method = "lm", se=FALSE, col='darkgrey') + labs(color="Merian element") +
  annotate("text", x = 10, y = 53, label = paste("R =", round(length_vs_repeat_corr, 2)), size = 3, hjust = 0) +
  annotate("text", x = 10, y = 52, label = paste("p =", round(length_vs_repeat_pval, 6)), size = 3, hjust = 0) +
  theme(legend.position = "none")


# plot length in P. icarus vs number chr in P. atlantica
Patlantica_Merians_accounting_fusions <- Patlantica_chr_assignments
Patlantica_Merians_accounting_fusions$Merian_element <- gsub(",", "+", Patlantica_Merians_accounting_fusions$Merian_element)

Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M17', 'M20'), 'M17+M20', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M25', 'M5'), 'M25+M5', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M14', 'M26'), 'M14+M26', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M18', 'M30'), 'M18+M30', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M1', 'M19'), 'M1+M19', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M28', 'M31'), 'M28+M31', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M12', 'M29'), 'M12+M29', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('M11', 'M23'), 'M11+M23', Patlantica_Merians_accounting_fusions$Merian_element)
Patlantica_Merians_accounting_fusions$Merian_element <- ifelse(Patlantica_Merians_accounting_fusions$Merian_element %in% c('MZ', 'M24'), 'M24+MZ', Patlantica_Merians_accounting_fusions$Merian_element)

# merge with P. icarus assignments that haven't undergone transformation to long format
temp <- merge(Picarus_chr_assignments, Picarus_lengths, by.x="query_chr",by.y="scaff")[,c(1,3,4)]
colnames(temp) <- c('scaff_Picarus', 'Merian_element', 'length_Picarus')
temp$Merian_element <- gsub(",", "+", temp$Merian_element)
Patlantica_Merians_accounting_fusions <- merge(Patlantica_Merians_accounting_fusions, temp, by='Merian_element')


# re-calculate number of chr per merian
counts2lengths <- Patlantica_Merians_accounting_fusions %>%
  group_by(Merian_element) %>%
  summarise(count = n())
lengths <- Patlantica_Merians_accounting_fusions[,c(1,5)] %>% unique()

counts2lengths <- merge(counts2lengths,lengths, by="Merian_element" )
counts2lengths$length_Picarus <- counts2lengths$length_Picarus / 1000000

counts2lengths$chr_type <- "Autosome"
counts2lengths$chr_type[counts2lengths$Merian_element == 'M24+MZ'] <- 'Z1'
counts2lengths$chr_type[counts2lengths$Merian_element == 'M16'] <- 'Z2'
counts2lengths_autosomes <- counts2lengths[counts2lengths$chr_type == "Autosome",]
counts2lengths_Z1 <- counts2lengths[counts2lengths$Merian_element == "M24+MZ",]
counts2lengths_Z2 <- counts2lengths[counts2lengths$Merian_element == "M16",]


# check for correlation between number of chr in P. atlantica and length of the chr in P. icarus
length_vs_count_corr_test <- cor.test(counts2lengths_autosomes$length_Picarus, counts2lengths_autosomes$count, method = "spearman")
length_vs_count_corr <- length_vs_count_corr_test$estimate # R: -0.2651694
length_vs_count_pval <- length_vs_count_corr_test$p.value #p-val 3.834381e-05
cat('R:', length_vs_count_corr, 'P-val:', length_vs_count_pval, '\n')
# conclude: strong correlation!

# plot repeat density vs chr length in P. icarus
repeats_vs_length_in_Picarus <- 
  ggplot(repeats_df_autosomes, aes(x=length_Picarus, y=repeat_per)) + 
  geom_point(size = 2, aes(color=factor(Merian_element_f))) +
  labs(x="Chr length in P. icarus (Mb)", y = "Repeat density in P. atlantica (%)") + 
  theme_classic() +  ggtitle('B') + theme(plot.title = element_text(face = "bold")) +
  geom_smooth(method = "lm", se=FALSE, col='darkgrey') + labs(color="Merian element") +
  annotate("text", x = 10, y = 82, label = paste("R =", round(length_vs_repeat_corr, 2)), size = 3, hjust = 0) +
  annotate("text", x = 10, y = 81, label = paste("p =", round(length_vs_repeat_pval, 6)), size = 3, hjust = 0) +
  theme(legend.background = element_rect(fill = "white", color = "darkgrey")) 

# + theme(legend.position = "none")

repeats_vs_length_in_Picarus
# plot repeat density vs chr length in P. atlantica
repeats_vs_length_in_Patlantica <- 
  ggplot(repeats_df_autosomes, aes(x=length_mb, y=repeat_per, color=Merian_element_f)) + geom_point()  +
  labs(x="Chr length (Mb)", y = "Repeat density (%)") + labs(color="Merian element") +
  theme_classic()  + theme(plot.title = element_text(face = "bold")) + xlim(0,6)

# plot number of chr in P. atlantica vs length in P. icarus 
col_palette <- c('#3d405b', '#e07a5f', '#81b29a')
Patlantica_chr_vs_Picarus_length_plot <- ggplot(counts2lengths_autosomes, aes(x=length_Picarus, y=count)) +
  geom_point(size = 2, aes(color=factor(chr_type))) +
  scale_color_manual(values=col_palette)+
  labs(x="Chr length in P. icarus (Mb)", y = "Number of chromosomes in P. atlantica") + 
  theme_classic() +  ggtitle('A') + theme(plot.title = element_text(face = "bold")) +
  geom_smooth(method = "lm", se=FALSE, col='darkgrey') + labs(color="Chromosome type")  +
  geom_point(data=counts2lengths_Z1, aes(length_Picarus, y=count, color=factor(chr_type))) +
  geom_point(data=counts2lengths_Z2, aes(length_Picarus, y=count, color=factor(chr_type))) +
  annotate("text", x = 10, y = 17, label = paste("R =", round(length_vs_count_corr, 2)), size = 3, hjust = 0) +
  annotate("text", x = 10, y = 16.5, label = paste("p =", round(length_vs_count_pval, 6)), size = 3, hjust = 0) +
  theme(legend.position = c(0.85,0.3),
        legend.background = element_rect(fill = "white", color = "darkgrey"))
Patlantica_chr_vs_Picarus_length_plot
# Plot histogram of chr lengths (autosomes only)
repeats_df$chr_type <- 'Autosome'
repeats_df$chr_type[repeats_df$scaff %in% c('SUPER_Z1', 'SUPER_Z2', 'SUPER_W1', 'SUPER_W2')] <- 'Sex_chr'

# Plot histogram of gene counts (autosomes only)
genes_df$chr_type <- 'Autosome'
genes_df$chr_type[genes_df$scaff %in% c('SUPER_Z1', 'SUPER_Z2', 'SUPER_W1', 'SUPER_W2')] <- 'Sex_chr'
genes_df_autosomes <- genes_df[genes_df$chr_type == 'Autosome',]

cat('Average number of genes per autosome:', mean(genes_df_autosomes$gene_count))
cat('SD in the number of genes per autosome:', sd(genes_df_autosomes$gene_count))


histogram_gene_counts <- ggplot(genes_df_autosomes, aes(x=gene_count)) + geom_histogram() + theme_classic() + 
  labs(x="Gene count", y="Number of chromosomes")  + geom_histogram(fill="lightgrey", color="black") +
  geom_vline(aes(xintercept=mean(gene_count)), color="black", linetype="dashed", size=1) + 
  ggtitle('C') + theme(plot.title = element_text(face = "bold"))

histogram_gene_counts

histogram_autosome_lengths <- ggplot(repeats_df_autosomes, aes(x=length_mb)) + geom_histogram() + theme_classic() + 
  labs(x="Chromosome length (Mb)", y="Number of chromosomes") + xlim(0,6) + geom_histogram(fill="lightgrey", color="black") +
  geom_vline(aes(xintercept=mean(length_mb)), color="black", linetype="dashed", size=1) + 
  ggtitle('B') + theme(plot.title = element_text(face = "bold"))

# plot autosomes vs sex chr - axis hugely distorted by sex chr..
ggplot(repeats_df, aes(x=length_mb, color=chr_type)) +  theme_classic() + 
  labs(x="Chr length in P. atlantica (Mb)", y="Number of chromosomes") + xlim(0,70) +
  geom_histogram(fill="white", alpha=0.5, position="identity", binwidth=0.5) 

# get mean autosome length:
mean(repeats_autosomes$length_mb)
sd(repeats_autosomes$length_mb)
# save output
fig1_B <- histogram_autosome_lengths
fig1_C <- histogram_gene_counts
fig1_B_C <- fig1_B + fig1_C
fig1_B_C
sup_fig_2 <- repeats_vs_length_in_Patlantica
fig2_A_B <- Patlantica_chr_vs_Picarus_length_plot + repeats_vs_length_in_Picarus 
fig2_A_B

#ggsave(fig2_A_B, filename='merian_legend.pdf', width=2, height=2)
ggsave(sup_fig_2, filename='repeat_density_vs_chr_length_Patlantica.pdf', device = 'pdf', width = 8, height = 5)
ggsave(sup_fig_2, filename='repeat_density_vs_chr_length_Patlantica.png', device = 'png', width = 8, height = 5)

ggsave(fig2_A_B, filename='chr_length_Picarus_vs_repeat_density_Merian_counts_220824.pdf', device = 'pdf', width = 12, height = 5)
ggsave(fig2_A_B, filename='chr_length_Picarus_vs_repeat_density_Merian_counts_220824.png', device = 'png', width = 12, height = 5.2)

ggsave(fig1_C, filename='histogram_autosome_lengths.pdf', device = 'pdf', width = 6, height = 5)
ggsave(fig1_C, filename='histogram_autosome_lengths.png', device = 'png', width = 6, height = 5.2)

ggsave(fig1_D, filename='histogram_gene_counts.pdf', device = 'pdf', width = 6, height = 5)
ggsave(fig1_D, filename='histogram_gene_counts.png', device = 'png', width = 6, height = 5.2)

ggsave(fig1_B_C, filename='histogram_chr_lengths_and_gene_counts.pdf', device = 'pdf', width = 7, height = 4)
ggsave(fig1_B_C, filename='histogram_chr_lengths_and_gene_counts.png', device = 'png', width = 7, height = 4)

