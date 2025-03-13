library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsignif)   
library(ggpubr)


# Aim: plot ds values per segment and across chr
setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/2_sex_chr")

locs <- read.csv('BRAKER_1_BRAKER_2_TSEBRA.output.transcript_locations.tsv', sep=' ', header=FALSE)
W1_Z1_dS <- read.csv('W1_Z1.codeml_NG_dS_results.reformatted.txt', sep=' ', header=FALSE)
W1_Z2_dS <- read.csv('W1_Z2.codeml_NG_dS_results.reformatted.txt', sep=' ', header=FALSE)
W2_Z2_dS <- read.csv('W2_Z2.codeml_NG_dS_results.reformatted.txt', sep=' ', header=FALSE)
best_hits <- read.csv('summarised_repicocal_blast_hits.txt', sep=' ', header=FALSE)

colnames(best_hits) <- c('rest_protein', 'rest_chr', 'W_protein', 'W_chr')
colnames(W1_Z1_dS) <- c('protein', 'dS')
colnames(W1_Z2_dS) <- c('protein', 'dS')
colnames(W2_Z2_dS) <- c('protein', 'dS')

colnames(locs) <- c('protein', 'chr', 'start', 'stop')
locs$midpos <- abs(locs$stop + locs$start)/ 2
locs$midpos_mb <- locs$midpos / 1000000
head(locs)

W1_Z1_dS <- merge(W1_Z1_dS, locs, by='protein')
W1_Z2_dS <- merge(W1_Z2_dS, locs, by='protein')
W2_Z2_dS <- merge(W2_Z2_dS, locs, by='protein')
head(W1_Z1_dS)

# filter clear outlier from W2_Z2
W1_Z2_dS <- W1_Z2_dS[W1_Z2_dS$dS > -0.2,] # removes the one outlier with ds=-1.
W1_Z2_dS <- W1_Z2_dS[W1_Z2_dS$dS <= 0.4,] # removes the one outlier with ds=0.55

# lets get average dS per chr pair
print(paste('Average dS W1 Z1:', mean(W1_Z1_dS$dS)*100)) # 7.2%
print(paste('Average dS W1 Z2:', mean(W1_Z2_dS$dS)*100)) # 9.4%
print(paste('Average dS W2 Z2:', mean(W2_Z2_dS$dS)*100)) # 1.8%

# plot per segment
alpha_val = 0.2 # set transparency level
W1_Z1_plot <- ggplot(data=W1_Z1_dS, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val) + theme_bw() + ggtitle('W1_Z1') +
  labs(x='Protein pos on W1 (Mb)', y='Synonymous site divergence (%)') +
  ylim(-0.1, 40)

W1_Z2_plot <- ggplot(data=W1_Z2_dS, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val) + theme_bw()+ ggtitle('W1_Z2') +
  labs(x='Protein pos on W1 (Mb)', y='Synonymous site divergence (%)') +
  ylim(-0.1, 40)

W2_Z2_plot <- ggplot(data=W2_Z2_dS, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val) + theme_bw() + ggtitle('W2_Z2') +
  labs(x='Protein pos on W2 (Mb)', y='Synonymous site divergence (%)') +
  ylim(-0.1, 40)

W1_Z1_plot + W1_Z2_plot + W2_Z2_plot + plot_layout(guides='collect')


# lets plot the two sets of proteins on W1
W1_Z2_dS$Z_chr <- 'Z2a'
W1_Z1_dS$Z_chr <- 'Z1'
W2_Z2_dS$Z_chr <- 'Z2b'
W1_vs_Zs <- rbind(W1_Z1_dS, W1_Z2_dS)

alpha_val =0.6
colour_palette <- c('#fe6d73','#ffcb77','#17c3b2')
W1_vs_Zs_plot <-ggplot(data=W1_vs_Zs, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val, aes(colour = Z_chr)) + theme_bw() + ggtitle('W1') +
  labs(x='Chr position (Mb)', y='Synonymous site divergence (%)') +
  scale_color_manual(values=colour_palette) +  
  theme(plot.title = element_text(face = "bold")) +
  ylim(-0.1, 55)

colour_palette <- c('#17c3b2')
W2_Z2_plot <- ggplot(data=W2_Z2_dS, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val, aes(colour = Z_chr))  + theme_bw() + ggtitle('W2') +
  labs(x='Chr position (Mb)', y='Synonymous site divergence (%)') +
  scale_color_manual(values=colour_palette) +
  theme(plot.title = element_text(face = "bold")) +
  ylim(-0.1, 55) + theme(axis.title.y = element_blank())

W1_vs_Zs_plot + W2_Z2_plot + plot_layout(guides='collect')

# lets make a boxplot of the dS values per chr pair
W1_Z1_dS$chr_pair <- 'W1_Z1'
W1_Z2_dS$chr_pair <- 'W1_Z2'
W2_Z2_dS$chr_pair <- 'W2_Z2'

all_dS <- rbind(W1_Z1_dS, W1_Z2_dS, W2_Z2_dS)

# use same colours as dominant Merian for consistency within figure
col_palette <-c('#fe6d73','#ffcb77', '#17c3b2')
box_plot_per_combo <- ggplot(all_dS, aes(x=chr_pair, y=dS*100)) + 
  geom_boxplot(aes(fill=chr_pair)) + theme_bw() +
#  geom_jitter(color="darkgrey", size=0.4, alpha=0.5) + # uncomment to see raw data points
  labs(x='Chr pair', y='Synonymous site divergence (%)') + 
  scale_fill_manual(values=col_palette) +
  theme(legend.position = "none") + 
  geom_signif(comparisons=list(c("W1_Z1", "W1_Z2"), 
                               c("W1_Z1", "W2_Z2")),
              map_signif_level = TRUE, y_position=c(38, 40))
box_plot_per_combo
# geom_signif uses a wilcox.test

col_palette <-c('#fe6d73','#ffcb77', '#17c3b2')
histogram_per_combo <- ggplot(all_dS, aes(x=dS*100, fill=chr_pair), ) +
  geom_histogram(color='black', alpha=0.6, bins=30) +
  scale_fill_manual(values=col_palette) + theme_bw() + 
  facet_grid(. ~ chr_pair) + coord_flip() + labs(y='Count', x='Synonymous site divergence (%)')
histogram_per_combo

ggsave(box_plot_per_combo, filename='boxplot_synonymous_site_divergence_per_sex_chr_pair.pdf', height=8, width=8)
ggsave(histogram_per_combo, filename='histogram_of_dS_per_chr_pair.pdf')

# Plot positions of proteins along chr and histograms together
panel_C_plus_D <- W1_vs_Zs_plot + W2_Z2_plot + plot_layout(guides='collect') + histogram_per_combo
panel_C_plus_D

ggsave(panel_C_plus_D, filename='dS_per_pair_along_chr_histograms_of_dS_plot.pdf', width=13, height=5)

cat(paste('Num of orthologues plotted:', nrow(W1_Z1_dS) + nrow(W1_Z2_dS) + nrow(W2_Z2_dS)))

      
# lets really look at where the super lowly diverged genes are:
W1_Z1_dS$high_or_low <- 'high'
W1_Z1_dS$high_or_low[W1_Z1_dS$dS < 0.02] <- 'low'

W1_Z1_dS_filt <- W1_Z1_dS[W1_Z1_dS$dS <= 0.02,]
W1_Z1_dS_filt_plot <- ggplot(data=W1_Z1_dS_filt, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val) + theme_bw() + ggtitle('W1_Z1') +
  labs(x='Protein pos on W1 (Mb)', y='Synonymous site divergence (%)') 
W1_Z1_dS_filt_plot


# lets look along Z chr instead of W chr

temp <- W1_Z1_dS[,c(1,2)]
colnames(temp) <- c('W_protein', 'dS')
W1_Z1_dS_Zlocs <- merge(temp, best_hits)
head(W1_Z1_dS_Zlocs)
W1_Z1_dS_Zlocs <- W1_Z1_dS_Zlocs[,c(3,4,2)]
colnames(W1_Z1_dS_Zlocs) <- c('protein', 'Z_chr', 'dS')
W1_Z1_dS_Zlocs <- merge(W1_Z1_dS_Zlocs, locs, by='protein')
head(W1_Z1_dS_Zlocs)

W1_Z1_dS_Zlocs$high_or_low <- 'high'
W1_Z1_dS_Zlocs$high_or_low[W1_Z1_dS_Zlocs$dS < 0.02] <- 'low'
Z_high_low_plot <- ggplot(data=W1_Z1_dS_Zlocs, color=high_or_low, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val, aes(color=high_or_low)) + theme_bw() + ggtitle('Z1') +
  labs(x='Protein pos on Z1 (Mb)', y='Synonymous site divergence (%)') 
Z_high_low_plot
W_high_low_plot <- ggplot(data=W1_Z1_dS, color=high_or_low, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val, aes(color=high_or_low)) + theme_bw() + ggtitle('W1') +
  labs(x='Protein pos on W1 (Mb)', y='Synonymous site divergence (%)') 
Z_high_low_plot + W_high_low_plot + plot_layout(guides='collect')

Z1_end <- W1_Z1_dS_Zlocs[W1_Z1_dS_Zlocs$midpos_mb >= 42.6,]
Z1_end_low <- Z1_end[Z1_end$high_or_low == 'low ',]
nrow(Z1_end_low) /nrow(Z1_end) # 55% of loci in last 3.2 Mb of Z1 are <2% divergence

Z1_rest <- W1_Z1_dS_Zlocs[W1_Z1_dS_Zlocs$midpos_mb < 42.6,]
Z1_rest_low <- Z1_rest[Z1_rest$high_or_low == 'low',]
nrow(Z1_rest_low) /nrow(Z1_rest) # 10% of loci on rest of chr have <2% divergence

W1_Z1_dS_Zlocs_filtered_last_3.2Mb <- W1_Z1_dS_Zlocs[W1_Z1_dS_Zlocs$midpos_mb < 42.6,]
mean(W1_Z1_dS_Zlocs_filtered_last_3.2Mb$dS)
mean(W1_Z2_dS$dS)

W1_Z1_dS_just_first_5Mb_W1 <- W1_Z1_dS[W1_Z1_dS$midpos_mb < 5,]
W1_Z1_dS_just_first_5Mb_W1_low <- W1_Z1_dS_just_first_5Mb_W1[W1_Z1_dS_just_first_5Mb_W1$high_or_low == 'low',]
nrow(W1_Z1_dS_just_first_5Mb_W1_low) / nrow(W1_Z1_dS_just_first_5Mb_W1) # 56.7%
W1_Z1_dS_rest <- W1_Z1_dS[W1_Z1_dS$midpos_mb > 5,]
W1_Z1_dS_rest_low <- W1_Z1_dS_rest[W1_Z1_dS_rest$high_or_low == 'low',]
nrow(W1_Z1_dS_rest_low) / nrow(W1_Z1_dS_rest) # 11%

W1_Z1_dS_without_first_5Mb_W1 <- W1_Z1_dS[W1_Z1_dS$midpos_mb > 5,]
mean(W1_Z1_dS_without_first_5Mb_W1$dS) # 0.096
mean(W1_Z1_dS$dS) # 0.072
mean(W1_Z2_dS$dS) # 0.094

# lets plot Z1 genes but colour by whether they are in the first/last 5 Mb of W1
W1_Z1_dS$Region[W1_Z1_dS$midpos_mb <= 5] <-  '< 5 Mb'
W1_Z1_dS$Region[W1_Z1_dS$midpos_mb > 5] <-  '> 5 Mb'
temp <- subset(W1_Z1_dS, select=c('protein', 'Region'))
temp <- merge(best_hits, temp, by.x='W_protein', 'protein')
temp <- subset(temp, select=c('rest_protein', 'Region'))
W1_Z1_dS_Zlocs <- merge(W1_Z1_dS_Zlocs, temp, by.x='protein', by.y='rest_protein')

# make plot
col_palette <- c('#fe6d73','#FFAEB1')
Z_chr_annotated_W1_dS <- ggplot(data=W1_Z1_dS_Zlocs, color=Region, aes(x=midpos_mb, y=dS*100)) + 
  geom_point(alpha=alpha_val, aes(color=Region)) + theme_bw() + ggtitle('Z1') +
  labs(x='Chr pos (Mb)', y='Synonymous site divergence (%)') +
  scale_color_manual(values=col_palette) +
  theme(legend.position = "inside") + 
  theme(legend.position.inside = c(0.2,0.85),
        legend.background = element_rect(fill="white", color="grey"))

# combine the two panels
Z_chr_annotated_W1_dS_and_all_Z_genes_plot <- Z_chr_annotated_W1_dS +
ggtitle('A') + theme(plot.title = element_text(face="bold")) + ggtitle('B') +
  theme(plot.title = element_text(face="bold"))

Z_chr_annotated_W1_dS_and_all_Z_genes_plot

ggsave(Z_chr_annotated_W1_dS_and_all_Z_genes_plot, filename='Z_chr_annotated_W1_dS_and_all_Z_genes_plot.png', width=10, height=5)
ggsave(Z_chr_annotated_W1_dS_and_all_Z_genes_plot, filename='Z_chr_annotated_W1_dS_and_all_Z_genes_plot.pdf', width=10, height=5)



#lets check distribution of all Z1-linked genes
gene_cov <- read.csv('../1_genomes/Polyommatus_atlantica.braker1_and_braker2.TSEBRA.gene_density_per_100kb.tsv', header=FALSE, sep='\t')[,c(1,2,3,7)]
colnames(gene_cov) <- c('chr', 'start', 'stop', 'prop')
gene_cov$midpos <- (gene_cov$start + gene_cov$stop)/2
head(gene_cov)
Z1_cov <- gene_cov[gene_cov$chr == 'SUPER_Z1',]

Z_genes_plot <- ggplot(data=Z1_cov, aes(x=midpos/1000000, y=prop*100)) + geom_point() + theme_bw() + labs(x='Chr position (Mb)',y='Gene density (%)')

c
ggsave(Z_genes_plot, filename='gene_density_per_100kb_on_Zchr.pdf', width=10, height=10)
ggsave(Z_genes_plot, filename='gene_density_per_100kb_on_Zchr.png', width=6, height=6)

length_diff <- read.csv('difference_in_coding_sequences_on_Zs_vs_Ws.tsv', sep='\t', header=FALSE)
colnames(length_diff) <- c('rest_protein', 'length_diff')
length_diff <- merge(length_diff, best_hits, by='rest_protein')
length_diff <- merge(length_diff, all_dS, by.x='W_protein', by.y='protein')
length_diff$Region <- '> 5 Mb'
length_diff$Region[length_diff$midpos_mb <= 5] <-  '< 5 Mb'
length_diff$Region[length_diff$chr_pair =="W2_Z2"] <-  'NA'
length_diff$Region[length_diff$chr_pair =="W1_Z2"] <-  'NA'

col_palette <- c('#fe6d73', '#FFAEB1', 'darkgrey')
length_diff_plot <- ggplot(length_diff, aes(x=length_diff, fill=Region)) + 
  geom_histogram(bins=15) + fill_palette(col_palette) +
  facet_grid(. ~ chr_pair) + theme_bw() + labs(x="Difference in W vs Z coding sequence length (nt)",y='Count')

ggsave(length_diff_plot, filename='Difference_in_W_vs_Z_sequence_length.png', width=10, height=4)

P_atlantica_vs_P_icarus_dS <- read.csv('P_atlantica_vs_P_icarus_M24.codeml_NG_dS_results.reformatted.txt', sep=' ', header=FALSE)
colnames(P_atlantica_vs_P_icarus_dS) <- c('P_icarus_M24', 'protein', 'dS_P_icarus')

P_atlantica_vs_P_icarus_dS <-merge(P_atlantica_vs_P_icarus_dS, W1_Z1_dS_Zlocs, by='protein')

# calculate the difference in dS between W1/Z1 in P.atlantica vs P.icarus vs P. atlantica
P_atlantica_vs_P_icarus_dS$diff_cf_Picarus <- P_atlantica_vs_P_icarus_dS$dS - P_atlantica_vs_P_icarus_dS$dS_P_icarus 

# lets add the label of whether the protein is on the first 5 Mb of W1 or not
temp <- subset(length_diff, select = c('rest_protein', 'Region'))
colnames(temp) <- c('protein', 'region_of_W')
P_atlantica_vs_P_icarus_dS <- merge(P_atlantica_vs_P_icarus_dS, temp, by='protein')

col_palette <- c('#fe6d73','#FFAEB1')
difference_in_dS_between_chr_vs_species <- ggplot(data=P_atlantica_vs_P_icarus_dS, aes(x=midpos_mb, y=diff_cf_Picarus*100, color=region_of_W)) +
  geom_point() + theme_bw() + labs(x='Chr position (Mb)', y='dS(P. atlantica Z1 & W1) - (dS(P. atlantica & P. icarus) (%)') + geom_smooth(se = FALSE) +
  scale_color_manual(values=col_palette)
difference_in_dS_between_chr_vs_species

ggplot(data=P_atlantica_vs_P_icarus_dS, aes(x=midpos_mb, y=diff_cf_Picarus*100)) +
  geom_point(aes(color=region_of_W)) + theme_bw() + labs(x='Chr position (Mb)', y='(dS(P. atlantica & P. icarus) - dS(P. atlantica Z1 & W1) (%)') + geom_smooth()

col_palette <- c('#fe6d73','#FFAEB1')
boxplot_diff_P_icarus <- ggplot(P_atlantica_vs_P_icarus_dS, aes(x=region_of_W, y=diff_cf_Picarus*100)) + 
  geom_boxplot(aes(fill=region_of_W)) + theme_bw() +
  geom_jitter(color="black", size=0.4, alpha=0.5) + #  see raw data points
  labs(x='Region of W1', y='dS(P. atlantica & P. icarus) - dS(P. atlantica Z1 & W1) (%)') + 
  scale_fill_manual(values=col_palette) +
  theme(legend.position = "none") + 
  geom_signif(comparisons=list(c("< 5 Mb", "> 5 Mb")),
              map_signif_level = TRUE)

# plot positions and boxplots together
normalised_linear_and_boxplots <- difference_in_dS_between_chr_vs_species + ggtitle('A') + theme(plot.title = element_text(face="bold")) +
  boxplot_diff_P_icarus + ggtitle('B') +
  theme(plot.title = element_text(face="bold"))

ggsave(normalised_linear_and_boxplots, filename='difference_in_dS_between_W1Z1_vs_P_icarus_linear_and_histogram.png', width=10, height=5)
ggsave(normalised_linear_and_boxplots, filename='difference_in_dS_between_W1Z1_vs_P_icarus_linear_and_histogram.pdf', width=10, height=5)


ggsave(difference_in_dS_between_chr_vs_species, filename='difference_in_dS_between_W1Z1_vs_P_icarus.png', width=6, height=5)
ggsave(difference_in_dS_between_chr_vs_species, filename='difference_in_dS_between_W1Z1_vs_P_icarus.pdf', width=6, height=5)

ggsave(boxplot_diff_P_icarus, filename='boxplot_difference_in_dS_between_W1Z1_vs_P_icarus.png', width=6, height=5)
ggsave(boxplot_diff_P_icarus, filename='boxplot_difference_in_dS_between_W1Z1_vs_P_icarus.pdf', width=6, height=5)

# combine to make one supfig


layout <- "
AB
CD
EE
"

sup_fig_divergence_analysis_combined <- Z_chr_annotated_W1_dS + ggtitle('A') +  theme(plot.title=element_text(face='bold')) +
  Z_genes_plot + ggtitle('B') +  theme(plot.title=element_text(face='bold')) +
  difference_in_dS_between_chr_vs_species + ggtitle('C') + theme(plot.title=element_text(face='bold'), 
                                                                 legend.position = "none") +
  boxplot_diff_P_icarus + ggtitle('D')  +  theme(plot.title=element_text(face='bold')) +
  length_diff_plot + ggtitle('E') + theme(plot.title=element_text(face='bold')) +
  plot_layout(design=layout)

ggsave(sup_fig_divergence_analysis_combined, filename='sup_fig_divergence_analysis_combined.pdf', width=12, height=12)


