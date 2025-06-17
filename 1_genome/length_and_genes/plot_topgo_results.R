# this script is adapted from: https://github.com/lstevens17/heligmosomoides_MS/blob/main/4_divergent_haplotypes/go_term_enrichment/plot_go_terms.R

library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

csv <- read.csv(args[1], header=T,sep='\t')

p <- ggplot(csv, aes(x=reorder(Term, -classicFisher), y=-log2(classicFisher), size=Significant/Expected)) + 
  geom_point(fill="#0092c1", pch=21) + coord_flip() + 
  geom_hline(yintercept=-log2(0.05), linetype=2) + 
  ylab("-log2(p-value)") + xlab("GO term") + 
  theme_bw() + scale_size(breaks = c(2, 4, 6, 8), name="Fold enrichment") + 
  theme(legend.position = c(0.85, 0.2), legend.box.background = element_rect(color="black", size=0.5))

output_file <- paste0('goterm_enrichment.',args[2],'.png')
p
ggsave(output_file, plot = p, width=10, height=5, units="in")
