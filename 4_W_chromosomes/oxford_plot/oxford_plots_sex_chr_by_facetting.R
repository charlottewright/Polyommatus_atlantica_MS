

library(ggplot2)
library(tidyverse)
library(patchwork)
library(forcats)
library(ggforce)
setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/2_sex_chr')

read_duplicated_busco <- function(buscoFile){
  busco <- read_tsv(buscoFile,
                    col_names = c("busco_id", "Status", "Sequence",
                                  "start", "end", "strand", "Score", "Length",
                                  "OrthoDB_url", "Description"),
                    col_types = c("ccciicdicc"),
                    comment = "#") %>%
    filter(Status == "Duplicated") %>%
    select(busco_id, Sequence, start, end)
  
  return(busco)
}


make_pairwise_dataframe <- function(df1, df2){
  df1$midpos_x <- df1$midpos
  df2$midpos_y <- df2$midpos
  pairwise_df = merge(df1, df2, by="busco_id")
  return(pairwise_df)
}

filter_buscos_x <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.x) %>%
    mutate(nGenes = n(),
           mxGpos = max(midpos_x)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}

filter_buscos_y <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.y) %>%
    mutate(nGenes = n(),
           mxGpos_y = max(midpos_y)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}

source('../functions_oxford_plots.R')
Merian_assignments_ref <- read_busco('Merian_elements_full_table.tsv')
Merian_assignments_ref <- select(Merian_assignments_ref, -c(start, end))
colnames(Merian_assignments_ref) <- c('busco_id', 'Merian')

species_1 <- read_duplicated_busco('../1_genomes/Polyommatus_atlantica.tsv') # species on Y-axis
species_1$midpos <- (species_1$start + species_1$end)/2

W_chr <- species_1[species_1$Sequence %in% c('SUPER_W1', 'SUPER_W2'),]
Z_chr <- species_1[species_1$Sequence %in% c('SUPER_Z1', 'SUPER_Z2'),]


colnames(Z_chr) <- c('busco_id', 'chr_1', 'start_1', 'end_1')
colnames(W_chr) <- c('busco_id', 'chr_2', 'start_2', 'end_2')

Z_and_W_chr <- merge(Z_chr, W_chr, by='busco_id')

Z_and_W_chr <- merge(Z_and_W_chr, Merian_assignments_ref, by="busco_id")
Z1_W1_end <-  data.frame(busco_id = 'placeholder1', chr_1 = 'SUPER_Z1', start_1 =46229683, end_1=46229683, chr_2 ='SUPER_W1', start_2 =65809200, end_2=65809200, Merian='NA')
Z2_W2_end <-  data.frame(busco_id = 'placeholder2', chr_1 = 'SUPER_Z2', start_1 =25546431, end_1=25546431, chr_2 ='SUPER_W2', start_2 =4527220, end_2=4527220, Merian='NA')

Z_and_W_chr <- rbind(Z_and_W_chr, Z1_W1_end)
Z_and_W_chr <- rbind(Z_and_W_chr, Z2_W2_end)

cols <- c("MZ" = "#227C9D", "M21" = "#17c3b2",
        "M16" = "#ffcb77", "M24" = "#fe6d73",
        "M7" = "grey", 'NA'='white')


p <- Z_and_W_chr %>% 
  ggplot(., aes(x=start_2/1e6, y=start_1/1e6, color=Merian)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(chr_1~chr_2, scales="free",space="free") + theme_bw() + 
  scale_colour_manual(values=cols) + xlab("Chr position (Mb)") + 
  ylab("Chr position (Mb)")  +
  # Force each facet to start from 0 on both x and y axes
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank())



p



ggsave("Ws_vs_Zs_oxford_plot.241024.pdf", plot = p, width = 20, height = 18, units = "cm")
ggsave("Ws_vs_Zs_oxford_plot.241024.png", plot = p, width = 20, height = 18, units = "cm")
