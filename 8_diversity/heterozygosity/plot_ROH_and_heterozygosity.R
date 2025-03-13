
setwd('/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/6_heterozygosity')

#het_f = sorted(list(Path(args.g).rglob('*.txt')))
#list_superfamilies = get_species_superfamily(args.f)
#print(list_superfamilies)
#df = calculate_average_heterozygosity(het_f, list_superfamilies, args.mcov, args.w)
#plot_species_heterozygosity(df, args.o)
#print(df)

file <- 'heterozygosity_output/Lysandra_coridon.deepvariant.PASS.biallelic.repeat_filtered.q15.vaf02.vaf08.dp10.dp93.heterozygosity.txt'
mcov = 6000 # default
window = 10000 # default
species_order <- c('Polyommatus_atlantica','Polyommatus_icarus','Polyommatus_iphigenia',
                   'Lysandra_bellargus','Lysandra_coridon', 'Plebejus_argus',
                   'Aricia_agestis','Aricia_artaxerxes','Cyaniris_semiargus', 'Celastrina_argiolus', 'Glaucopsyche_alexis', 'Phengaris_arion', 'Lycaena_phlaeas')

het_df <- data.frame()

files <- list.files(path = 'heterozygosity_output/')

for (i in files){
  df <- read.csv(paste0('heterozygosity_output/',i), sep='\t', header=FALSE)
  species <-  strsplit(i, "\\.")[[1]][1]
  print(species)
  colnames(df) <-  c('chrom', 'start', 'end', 'ncov', 'nhet', 'snp_count')
  df <- df[df$ncov > mcov,]
  df$snp_count_scaled_window <- df$snp_count / window
  avg_sem_het <- df$snp_count_scaled_window
  snp_count_avg <- mean(avg_sem_het)
  new_row <- data.frame(Species=species, Heterozygosity=snp_count_avg)
  het_df <- rbind(het_df, new_row)
}
het_df$Species <- factor(het_df$Species, levels = species_order)


# Now plot ROHs
genome_size_file <- 'species_genome_info.txt'
genome_size <- read.csv(genome_size_file, sep='\t', header=FALSE)
colnames(genome_size) <- c('species', 'genome_size')
#roh_files = list(Path(args.d).rglob('*.insideROH.txt'))
roh_df <- data.frame()

file <- 'runs_of_homozygosity_output/Aricia_agestis.deepvariant.PASS.biallelic.repeat_filtered.q15.vaf02.vaf08.dp10.dp102noZW.insideROH.txt'
files <- list.files(path = 'runs_of_homozygosity_output/insideROH/')

for (i in files){
  df <- read.csv(paste0('runs_of_homozygosity_output/insideROH/',i), sep='\t', header=FALSE)
  colnames(df) <- c('chromosome', 'start', 'end', 'size', 'avg_snp_count')
  df$length <- (df$end - df$start)/ 1000000
  species <-  strsplit(i, "\\.")[[1]][1]
  #myroh[species_name, genus, genome].append(length)
  total_roh_len <- sum(df$length)
  total_roh_num <- nrow(df)
  # Based on their length, classify ROHs into short, medium, and long
  short_ROH_df <- df[(df$length <= 0.1),]
  medium_ROH_df <- df[(df$length > 0.1 & df$length < 1.0),]
  long_ROH_df <- df[(df$length > 1.0),]
  short_ROH <- nrow(short_ROH_df)
  medium_ROH <- nrow(medium_ROH_df)
  long_ROH <- nrow(long_ROH_df)
  sum_short <- sum(short_ROH_df$length)
  sum_medium <- sum(medium_ROH_df$length)
  sum_long <- sum(long_ROH_df$length)
  print(species)
  #genome_size_for_species <- genome_size[genome_size$species == species,]
#  genome_size_for_species <- genome_size_for_species$genome_size[1]
#  print(genome_size_for_species)
  new_row <- data.frame(Species=species, Category='All', num_roh=total_roh_num, len_roh=total_roh_len)
  roh_df <- rbind(roh_df, new_row)
  new_row <- data.frame(Species=species, Category='Short', num_roh=short_ROH, len_roh=sum_short)
  roh_df <- rbind(roh_df, new_row)
  new_row <- data.frame(Species=species,  Category='Medium', num_roh=medium_ROH, len_roh=sum_medium)
  roh_df <- rbind(roh_df, new_row)
  new_row <- data.frame(Species=species,  Category='Large', num_roh=long_ROH, len_roh=sum_long)
  roh_df <- rbind(roh_df, new_row)
}

col_palette <- c('#FF715B','#E8C547', '#5C80BC')
roh_df$Species <- factor(roh_df$Species, levels = species_order)
roh_df <- roh_df[roh_df$Category != 'All',]

# Stacked barchart of number of ROHs per species, grouped by size
num_plot <- ggplot(roh_df, aes(fill=Category, y=num_roh, x=Species)) + theme_bw() +
  geom_bar(position="stack", stat="identity")  +
  labs(x='Species', y='Number of ROHs') + theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_manual(values=col_palette) + ggtitle('C') +  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position = 'none') +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text= element_text(size=10), axis.title=element_text(size=11)) 

# Stacked barchart of length of ROHs per species, grouped by size
length_plot <- ggplot(roh_df, aes(fill=Category, y=len_roh, x=Species)) + theme_bw() +
  geom_bar(position="stack", stat="identity") +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x='Species', y='Total length of ROHs') +   scale_fill_manual(values=col_palette) +
  ggtitle('B') +  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position = c(0.07, 0.65)) + theme(axis.title.x=element_blank(),
                                                 axis.text= element_text(size=10)) 


#length_plot_modified <-  length_plot + theme(axis.text.x= element_blank(),
#                               axis.ticks.x= element_blank(),
#                               axis.title.x=element_blank()) 

het_plot <- ggplot(data=het_df, aes(x=Species, y=Heterozygosity)) +geom_bar(stat = "identity", fill='#713E5A') + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ggtitle('A') +  theme(plot.title = element_text(face = "bold")) + theme(axis.title.x=element_blank(),
                                                                            axis.text= element_text(size=10)) 


het_plot_modified <- het_plot + theme(axis.title.x=element_blank()) 

combined_plot <- het_plot_modified + length_plot + num_plot + plot_layout(ncol=1) 
combined_plot

ggsave(combined_plot, filename='heteorzygosity_and_rohs_per_species.241024.pdf')
ggsave(combined_plot, filename='heteorzygosity_and_rohs_per_species.241024.png')
