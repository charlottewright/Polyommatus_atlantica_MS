library(ggplot2)
library(dplyr)

# Aim: Plot historgams of repeat density profiles of random vs breakpoint regions
print(random_file_name)
counter = 1
col_names = c('chr', 'repeat_density')

# Read in random regions
read_and_combine_random_files <- function(file_dir, file_prefix, repeat_types) {
  col_names = c('chr', 'repeat_density')
  counter <- 1
  result_df <- NULL
  for (i in repeat_types) {
    file_name <- paste0(file_prefix, i, '.gff.tsv')
    full_file_path <- file.path(file_dir, file_name)
    print(full_file_path)
    temp_df <- read.csv(full_file_path, sep='\t', header=FALSE)[, c(1, 7)]
    colnames(temp_df) <- col_names
    temp_df$repeat_type <- i
    if (counter == 1) {
      result_df <- temp_df
    } else {
      result_df <- rbind(temp_df, result_df)
    }
    counter <- counter + 1
  }
  result_df$repeat_density <- result_df$repeat_density*100 # convert to %
  return(result_df)
}
# Read in breakpoint regions
read_and_combine_breakpoint_files <- function(file_dir, file_prefix, repeat_types) {
  col_names = c('chr', 'repeat_density')
  counter <- 1
  result_df <- NULL
  for (i in repeat_types) {
    file_name <- paste0(file_prefix, i, '.gff.tsv')
    full_file_path <- file.path(file_dir, file_name)
    print(full_file_path)
    temp_df = read.csv(full_file_path, sep='\t', header=FALSE)[,c(1,7)]
    colnames(temp_df) <- col_names
    temp_df$repeat_type <- i
    if (counter == 1) {
      result_df <- temp_df
    } else {
      result_df <- rbind(temp_df, result_df)
    }
    counter <- counter + 1
  }
  result_df$repeat_density <- result_df$repeat_density*100 # convert to %
  return(result_df)
}
# Read in Lysandra coridon
setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/4_fissions")
repeat_types = c('LINE', 'SINE', 'DNA', 'RC', 'LTR', 'Repeats')

# Read in Polyommatus atlantica (currently using ilPolAtla version 3)
file_dir = 'temp_old_ilPolAtla/'
file_prefix <- 'repeat_density_of_random_regions_10000_Cyaniris_semiargus.filtered'
ilPolAtal_v3_random_df <- read_and_combine_random_files(file_dir, file_prefix, repeat_types)
file_prefix <- 'repeat_density_of_ilPolAtal1_v3_vs_ilCyaSemi_breakpoints_noZchr_040124.tsv_Cyaniris_semiargus.filtered'
ilPolAtal_v3_breakpoints_df <- read_and_combine_breakpoint_files(file_dir, file_prefix, repeat_types)

avg_ilPolAtal_v3_breakpoints_df <- ilPolAtal_v3_breakpoints_df %>% group_by(repeat_type) %>% summarize(mean_value = mean(repeat_density))

# Reorder the levels of the 'repeat_type' factor variable
ilPolAtal_v3_random_df$repeat_type <- factor(ilPolAtal_v3_random_df$repeat_type, levels = repeat_types)

avg_ilPolAtal_v3_breakpoints_df$repeat_type <- factor(avg_ilPolAtal_v3_breakpoints_df$repeat_type, levels = repeat_types)


# Calculate confidence interval 
mean_random <- mean(ilPolAtal_v3_random_df$repeat_density)
sd_random <- sd(ilPolAtal_v3_random_df$repeat_density)
conf_int <- mean_random + c(-1.96, 1.96) * sd_random/sqrt(length(ilPolAtal_v3_random_df)) 

# work out 95% confidence interval for each type of repeat
conf_int_table <- ilPolAtal_v3_random_df %>%
  group_by(repeat_type) %>%
  summarise(
    mean_repeat_density = mean(repeat_density, na.rm = TRUE),
    sd_repeat_density = sd(repeat_density, na.rm = TRUE),
    n = n(),  # Number of observations in each group
    conf_int_lower = mean_repeat_density - 1.96 * (sd_repeat_density / sqrt(n)),
    conf_int_upper = mean_repeat_density + 1.96 * (sd_repeat_density / sqrt(n)),
    sd_lower = mean_repeat_density - sd_repeat_density,
    sd_higher = mean_repeat_density + sd_repeat_density)

# decided not to use confidence intervals (CIs) because sample size is so high that CIs are tiny - not informative
ggplot(data=ilPolAtal_v3_random_df, aes(x=repeat_density)) + geom_histogram() +
  geom_vline(data = avg_ilPolAtal_v3_breakpoints_df, aes(xintercept = mean_value), linetype = "dotted", color = "red", size=1.2) + 
  theme_bw() + xlab('Repeat density (%)') + ylab('Count') +
  geom_ribbon(aes(ymin = 0, ymax = Inf, 
                  xmin = conf_int[1], 
                  xmax = conf_int[2]), 
              fill = "gray80", alpha = 0.5)


hist_plot <- ggplot() + 
  theme_bw() + 
  geom_histogram(data = ilPolAtal_v3_random_df, aes(x = repeat_density),fill='lightgrey',color='black',bins=50) +  # Histogram
  geom_rect(data=conf_int_table, aes(xmin=sd_lower, xmax=sd_higher, 
                                     ymax=0, ymin=Inf), fill="blue",alpha=0.2) + 
  geom_vline(data = avg_ilPolAtal_v3_breakpoints_df, aes(xintercept = mean_value), 
             color = "red", linetype = "dashed") +  # Vertical line for mean
  facet_wrap(~ repeat_type, scales="free") + labs(x='Repeat density (%)',y='Number of random regions')

ggsave(hist_plot, filename='repeat_density_random_vs_observed.ilPolAtlav3.pdf',height=8, width=12)
ggsave(hist_plot, filename='repeat_density_random_vs_observed.ilPolAtlav3.png',height=8, width=12)

  