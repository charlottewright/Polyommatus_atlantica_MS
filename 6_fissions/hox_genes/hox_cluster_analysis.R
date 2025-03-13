library(ggplot2)
library(gggenes)
library(dplyr)

setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/4_fissions")

HOX_cluster <- read.csv("HOX_cluster.reformatted.txt", header=FALSE, sep='')
PRD_cluster <- read.csv("PRD_cluster.reformatted.txt",header=FALSE, sep='')

colnames(HOX_cluster) <- c("chr", "species","start",	"end","gene","strand","percent_identity")
colnames(PRD_cluster) <- c("chr", "species","start",	"end","gene","strand","percent_identity")

HOX_cluster$species <- factor(HOX_cluster$species, levels = unique(HOX_cluster$species))
HOX_cluster$gene <- factor(HOX_cluster$gene, levels = unique(HOX_cluster$gene))

PRD_cluster$gene <- factor(PRD_cluster$gene, levels = unique(PRD_cluster$gene))

P_atlantica_HOX <- HOX_cluster[HOX_cluster$species== "Polyommatus_atlantica",]
P_icarus_HOX <- HOX_cluster[HOX_cluster$species== "Polyommatus_icarus",]
genes_of_interest <- c('lab', 'Abd-B', 'abd-A', 'Ubx', 'Antp', 'ftz', 'Scr', 'Dfd', 'zen', 'ShxB', 'ShxA', 'Pb', 'Ro', 'ShxC')

# just plot the core hox gene cluster
P_atlantica_HOX_filt <- P_atlantica_HOX %>% filter(gene %in% genes_of_interest)
P_icarus_HOX_filt <- P_icarus_HOX %>% filter(gene %in% genes_of_interest)
P_atlantica_HOX_filt$percent_identity <- as.integer(P_atlantica_HOX_filt$percent_identity)

# remove duplicate of ShxB which has lower identity
P_atlantica_HOX_filt <- P_atlantica_HOX_filt %>% filter(percent_identity > 55)
# remove duplicate of Ubx which has lower identity
P_atlantica_HOX_filt <- subset(P_atlantica_HOX_filt, !(gene == "Ubx" & percent_identity < 100))
# both Ro and Pb have two hits with 100% which are very close to eachother - due to being two exons of same gene (checked based on size of intron)
# so lets filter to just keep one exon per gene
P_atlantica_HOX_filt <- P_atlantica_HOX_filt %>%
  filter(!(gene == 'Ro' & start == 997098)) %>%
  filter(!(gene == 'Pb' & start == 1281769))

# remove duplicates of ShxC and zen which have lower identity
P_icarus_HOX_filt <- P_icarus_HOX_filt %>% filter(percent_identity > 55)
# remove duplicate of Ubx which has lower identity
P_icarus_HOX_filt <- subset(P_icarus_HOX_filt, !(gene == "Ubx" & percent_identity < 100))
# both Ro and Pb have two hits with 100% which are very close to eachother - due to being two exons of same gene (checked based on size of intron)
# so lets filter to just keep one exon per gene
P_icarus_HOX_filt <- P_icarus_HOX_filt %>%
  filter(!(gene == 'Ro' & start == 13403711)) %>%
  filter(!(gene == 'Pb' & start == 13162450))
 
P_atlantica_HOX_plot <- ggplot(P_atlantica_HOX_filt, aes(xmin =start/1000000, xmax = end/1000000, y = chr, fill = gene, label = gene, order= gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  facet_wrap(~ chr, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("lab" = "#65919b", "Abd-B" = "#dfc27d", "abd-A" = "#dfc27d", "Ubx" = "#dfc27d", "Antp" = "#dfc27d", "ftz" = "#dfc27d", "Scr" = "#dfc27d", "Dfd" = "#dfc27d", "zen" = "#dfc27d", "Pb" = "#dfc27d", "ShxD" = "#1f78b4", "ShxC" = "#33a02c", "ShxB" = "#e31a1c", "ShxA" = "#ff7f00", "Ro" = "#6a3d9a")) +
  theme_genes() +
  theme(axis.title.y=element_blank(), axis.text.y=element_text(size = 10)) + xlab('Chr position (Mb)')


P_icarus_HOX_plot <- ggplot(P_icarus_HOX_filt, aes(xmin = start/1000000, xmax = end/1000000, y = chr, fill = gene, label = gene, order= gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  facet_wrap(~ chr, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("lab" = "#65919b", "Abd-B" = "#dfc27d", "abd-A" = "#dfc27d", "Ubx" = "#dfc27d", "Antp" = "#dfc27d", "ftz" = "#dfc27d", "Scr" = "#dfc27d", "Dfd" = "#dfc27d", "zen" = "#dfc27d", "Pb" = "#dfc27d", "ShxD" = "#1f78b4", "ShxC" = "#33a02c", "ShxB" = "#e31a1c", "ShxA" = "#ff7f00", "Ro" = "#6a3d9a")) +
  theme_genes() +
  theme(axis.title.y=element_blank(), axis.text.y=element_text(size = 10)) + xlab('Chr position (Mb)')

P_atlantica_HOX_plot + P_icarus_HOX_plot + plot_layout(guides='collect')

ggsave(P_atlantica_HOX_plot, filename="Hox_genes_in_P_atlantica.png")
ggsave(P_atlantica_HOX_plot, filename="Hox_genes_in_P_atlantica.pdf")

ggsave(P_icarus_HOX_plot, filename="Hox_genes_in_P_icarus.png")
ggsave(P_icarus_HOX_plot, filename="Hox_genes_in_P_icarus.pdf")

write_tsv(P_icarus_HOX_filt, file ='P_icarus_HOX_gene_locations_filtered.tsv')
write_tsv(P_atlantica_HOX_filt, file ='P_atlantica_HOX_gene_locations_filtered.tsv')


# PRD analysis
P_atlantica_PRD <- PRD_cluster[PRD_cluster$species== "Polyommatus_atlantica",]
genes_of_interest <- c('Hbn', 'Rax', 'Arx', 'Otp')
P_atlantica_PRD_filt <- P_atlantica_PRD %>% filter(gene %in% genes_of_interest)
P_atlantica_PRD_filt <- P_atlantica_PRD_filt %>% filter(chr == "SUPER_83")

P_atlantica_PRD_plot <- ggplot(P_atlantica_PRD_filt, aes(xmin = (start+10000)/1000000, xmax = end/1000000, y = chr, fill = gene, label = gene, order= gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  facet_wrap(~ chr, scales = "free", ncol = 1) +
 # scale_fill_manual(values = c("lab" = "#65919b", "Abd-B" = "#dfc27d", "abd-A" = "#dfc27d", "Ubx" = "#dfc27d", "Antp" = "#dfc27d", "ftz" = "#dfc27d", "Scr" = "#dfc27d", "Dfd" = "#dfc27d", "zen" = "#dfc27d", "Pb" = "#dfc27d", "ShxD" = "#1f78b4", "ShxC" = "#33a02c", "ShxB" = "#e31a1c", "ShxA" = "#ff7f00", "Ro" = "#6a3d9a")) +
  theme_genes() + labs(x="Chr position (Mb)") + 
  theme(axis.title.y=element_blank(), axis.text.y=element_text(size = 10))


ggsave(P_atlantica_PRD_plot, filename="PRD_genes_in_P_atlantica.png")
ggsave(P_atlantica_PRD_plot, filename="PRD_genes_in_P_atlantica.pdf")

write_tsv(P_atlantica_PRD_filt, file ='P_atlantica_PRD_gene_locations_filtered.tsv')





