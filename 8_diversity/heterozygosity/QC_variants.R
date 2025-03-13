
library(dplyr)
library(ggplot2)
library(patchwork)

read_VCF <- function(species){
  vcf <- read.csv(paste('data/', species, '.deepvariant.PASS.biallelic.repeat_filtered.vcf.tsv', sep=''), sep='\t', skip=1, header=FALSE)[c(1,2,3,6,7)]
  colnames(vcf) <- c('chr', 'position', 'qual', 'DP', 'AD')
  vcf$VAF = vcf$AD / vcf$DP
  return(vcf)
}

filter_VCF <- function(vcf){
  vcf_filt <- vcf %>% group_by(chr) %>% mutate(max_pos = max(position)) %>% filter(max_pos > 600000) # > 0.5 Mb to remove schrapnel
  pre_filt <- unique(vcf$chr)
  post_filt <- unique(vcf_filt$chr)
  removed_scaffs <- setdiff(pre_filt, post_filt) # scaffolds removed
  print('Filtered out scaffolds:')
  print(removed_scaffs)
  print('Retained scaffolds: ')
  print(post_filt)
  return(vcf_filt)
}

setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/6_heterozygosity/")
make_QC_plots <- function(vcf_filt){
  # Look at DP - set max depth to be x2 that of the mean depth
  max_depth <- mean(vcf_filt$DP)*2
  # was 110
  A <- ggplot(data = vcf_filt, aes(x=DP, fill=chr)) + geom_density(alpha=0.4) + xlim(0,110) + theme_bw() +
    geom_vline(xintercept = 10, linetype="dotted") +
    geom_vline(xintercept = max_depth, linetype="dotted") +
    theme(legend.position = "none") #+ facet_wrap(~chr)# +
  
  # Look at QUAL - remove sites QUAL < 15
  # see discussion here https://github.com/freebayes/freebayes/issues/237
  B <- ggplot(data = vcf_filt, aes(x=qual, fill=chr)) + geom_density(alpha=0.4) + theme_bw() +xlim(0,100) + theme(legend.position = "none") +
    geom_vline(xintercept = 15, linetype="dotted") #+ facet_wrap(~chr)
  
  # Look at VAF - remove < 0.2 & > 0.8
  C <- ggplot(data = vcf_filt, aes(x=VAF, fill=chr)) + geom_density(alpha=0.4) + xlim(0,1) + theme_bw() +
    geom_vline(xintercept = 0.2, linetype="dotted") +
    geom_vline(xintercept = 0.8, linetype="dotted") +
    theme(legend.position = "none")
  
  plot_ABC <- A + B + C & theme(legend.position = "none") 
  plot_ABC <- plot_ABC + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")
  
  return(plot_ABC)
}

input_species <- "Aricia_agestis"  
input_species <- "Aricia_artaxerxes" 
input_species <- "Celastrina_argiolus"
input_species <- "Cyaniris_semiargus"
input_species <- "Glaucopsyche_alexis"
input_species <- "Lycaena_phlaeas"
input_species <- "Lysandra_bellargus"
input_species <- "Lysandra_coridon"
input_species <- "Plebejus_argus"
input_species <- "Polyommatus_atlantica"
input_species <- "Polyommatus_icarus"
 
species_list <- list('Aricia_agestis', 'Aricia_artaxerxes', 'Celastrina_argiolus', 'Cyaniris_semiargus', 
                     'Glaucopsyche_alexis', 'Lycaena_phlaeas', 'Lysandra_bellargus', 'Lysandra_coridon',
                     'Plebejus_argus', 'Polyommatus_atlantica', 'Polyommatus_icarus','Phengaris_arion')
# run code
vcf <- read_VCF(input_species)
vcf_filt <- filter_VCF(vcf)

#vcf_filt <- filter(vcf_filt, chr %in% c('SUPER_166', 'SUPER_167', 'SUPER_168', 'SUPER_201', 'SUPER_202', 'SUPER_203'))
species_list <- list('Polyommatus_iphigenia', 'Aricia_agestis', 'Aricia_artaxerxes')
n = 1
for (species in species_list) {
  print(species)
  vcf <- read_VCF(species)
  vcf_filt <- filter_VCF(vcf)
  plot_ABC <- make_QC_plots(vcf_filt)
  ggsave(plot=plot_ABC, filename=paste('figures/deepvariant_QC/', species, '_variant_QC.pdf', sep=''), width=12, height=5)
  ggsave(plot=plot_ABC, filename=paste('figures/deepvariant_QC/', species, '_variant_QC.png', sep=''), width=12, height=5)
}
#ggsave(plot=D, filename=paste('figures/deepvariant_QC/', input_species, '_variant_VAF_per_chr.png', sep=''), width=20, height=20)
#summary(vcf_filt$VAF)

input_species <- 'Polyommatus_iphigenia'
vcf <- read_VCF(input_species)
vcf_filt <- filter_VCF(vcf)
max_depth <- mean(vcf_filt$DP)*2
max_depth

