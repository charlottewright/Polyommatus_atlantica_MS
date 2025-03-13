
# Aim: compare simulations of fissioned chr to observed in P. atlantica
library(ggplot2)
library(patchwork)

setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/4_fissions/")
mean_sims <- read.csv('mean_BP_per_fragment_per_simulation_set.03032025.tsv', header=FALSE)
sd_sims <- read.csv('stdev_BP_per_fragment_per_simulation_set.03032025.tsv', header=FALSE)
colnames(mean_sims) <- 'mean_per_simulation'
colnames(sd_sims) <- 'stdev_per_simulation'
head(mean_sims)
head(sd_sims)

# mean for P. atlantica was 3.07
# stdev for P.atlantica was 1.73 
P_atlantica_mean = 3.07
P_atlantica_stdev = 1.73

mean_dist <- ggplot(mean_sims, aes(x=mean_per_simulation)) + geom_histogram(color="black", fill="white") + theme_bw() +
  geom_vline(xintercept = P_atlantica_mean, linetype = "dashed", color = "red") + # Add a horizontal line for P. atlantica mean 
  labs(x='Mean enriched BP per simulation',y='Count')

stdev_dist <- ggplot(sd_sims, aes(x=stdev_per_simulation)) + geom_histogram(color="black", fill="white") + theme_bw() +
  geom_vline(xintercept = P_atlantica_stdev, linetype = "dashed", color = "red") +  # Add a horizontal line for P. atlantics stdev
  labs(x='Stdev of enriched BP per simulation',y='Count')

combined_plot <- mean_dist + stdev_dist

mean(mean_sims$mean_per_simulation) # 3.2
mean(sd_sims$stdev_per_simulation) # 1.9


ggsave(combined_plot, filename='simulated_gene_enrichment_results.pdf', height=4, width=10)


