
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

# Aim: look at divergence landscape of repeat profile

setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/1_genomes")

read_in_div_file <- function(filename){
  df <- read.csv(filename, sep=' ', header=FALSE)
  colnames(df) <- c('div', 'chr', 'start', 'stop', 'type')
  df$Family <- sapply(str_split(df$type, "/"), `[`, 1)
  df$start <- df$start / 1000000
  df$stop <- df$stop / 1000000
  df$midpos <- (df$stop + df$start)/2
  head(df)
  df$Family[df$Family=='SINE?'] <- 'SINE'
  df$Family[df$Family=='Retroposon'] <- 'Unknown'
  return(df)
} 


col_palette <- c('#001219', '#005f73', '#0a9396', '#94d2bd', '#e9d8a6', '#ee9b00',
                 '#ca6702', '#DD4803', '#ae2012', '#751A1D',
                 'darkgrey','lightgrey')

Picarus <- read_in_div_file('Polyommatus_icarus.fa.earlGrey_divergence.tsv')
Patlantica <- read_in_div_file('Polyommatus_atlantica.fa.earlGrey_divergence.tsv')

#df <- read.csv('Polyommatus_atlantica.fa.earlGrey_divergence.tsv', sep=' ', header=FALSE)
#colnames(df) <- c('div', 'chr', 'start', 'stop', 'type')
#df$Family <- sapply(str_split(df$type, "/"), `[`, 1)
#df$start <- df$start / 1000000
#df$stop <- df$stop / 1000000
#head(df)
#df$Family[df$Family=='SINE?'] <- 'SINE'

# lets take a look at one chr for example
Patlantica_Z1 <- Patlantica[Patlantica$chr == 'SUPER_Z1',]
sex_chr <- Patlantica[Patlantica$chr %in% c('SUPER_Z1', 'SUPER_Z2', 'SUPER_W1', 'SUPER_W2'),]
autosomes <- Patlantica[!Patlantica$chr %in% c('SUPER_Z1', 'SUPER_Z2', 'SUPER_W1', 'SUPER_W2'),]

ggplot(Patlantica_Z1, aes(x=midpos, y=div, color=Family)) + geom_point() + theme_classic()

# look at just sex cgr
sex_chr_plot <- ggplot(sex_chr, aes(x=div, fill=Family)) + geom_histogram(col='black') + theme_classic() + 
  labs(x='Divergence (%)', y='Count') + scale_fill_manual(values=col_palette) + 
  ggtitle('D') + theme(plot.title = element_text(face = "bold")) + facet_grid(~chr)

# look at just autosomes
autosome_plot <- ggplot(autosomes, aes(x=div, fill=Family)) + geom_histogram(col='black') + theme_classic() + 
  labs(x='Divergence (%)', y='Count') + scale_fill_manual(values=col_palette) + 
  ggtitle('Autosomes') + theme(plot.title = element_text(face = "bold")) 

sex_chr_plot + autosome_plot +  plot_layout(guides = 'collect', nrow=2) 

Patlantica_repeat_div_plot <- ggplot(Patlantica, aes(x=div, fill=Family)) + geom_histogram(col='black') + theme_classic() + 
  labs(x='Divergence (%)', y='Count') + scale_fill_manual(values=col_palette) + 
  ggtitle('Polyommatus atlantica') + theme(plot.title = element_text(face = "bold"))

Picarus_repeat_div_plot <- ggplot(Picarus, aes(x=div, fill=Family)) + geom_histogram(col='black') + theme_classic() + 
  labs(x='Divergence (%)', y='Count') + scale_fill_manual(values=col_palette) + 
  ggtitle('Polyommatus icarus') + theme(plot.title = element_text(face = "bold"))

Patlantica_plus_Picarus_plot <- (Picarus_repeat_div_plot) / Patlantica_repeat_div_plot +  plot_layout(guides = 'collect')


ggsave(repeat_div_plot, filename='repeat_divergence_distribution.png', device = 'png', width = 6, height = 5.2)
ggsave(repeat_div_plot, filename='repeat_divergence_distribution.pdf', device = 'pdf', width = 6, height = 5.2)

ggsave(Patlantica_plus_Picarus_plot, filename='repeat_divergence_distribution_Patlantica_and_Picarus.png', device = 'png', width = 8, height = 10.4)
ggsave(Patlantica_plus_Picarus_plot, filename='repeat_divergence_distribution_Patlantica_and_Picarus.pdf', device = 'pdf', width = 8, height = 10.4)
