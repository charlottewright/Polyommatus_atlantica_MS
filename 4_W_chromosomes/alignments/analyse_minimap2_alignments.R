setwd("/Users/cw22/Documents/transfer_folder/GitHub_projects/polyommatus_atlantica/2_sex_chr")

# read in data

path_to_all_aln <-'W1_and_W2_vs_genome_asm20.paf'
path_to_all_aln_flipped <- 'genome_vs_W2_and_W2_asm20.paf'

#ali <- read_paf(path_to_aln)
all_ali <- read_paf(path_to_all_aln)
all_aln_flipped <- read_paf(path_to_all_aln_flipped)

# de = gap-compressed per-base sequence divergence
ggplot(all_ali, aes(alen, de)) + 
  geom_point(alpha=0.6, colour="steelblue", size=2) + 
  scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
  scale_y_continuous("Per base divergence") + 
  theme_pubr()

by_q <- aggregate(de ~ qname, data=all_ali, FUN=mean)
knitr::kable(by_q)

# don't need to filter secondary alignments as using option "secondary=none' in minimap2
nrow(all_ali) - nrow(prim_alignment)

long_ali_5000 <- subset(all_ali, alen > 5000 & mapq >= 40 )
all_long_ali_5000 <- subset(all_ali, alen > 5000 & mapq >= 40 )
all_flipped_long_ali_5000 <- subset(all_aln_flipped, alen > 5000 & mapq >= 40 )
nrow(long_ali_5000)

minimap2_aln <- dotplot(all_ali, label_seqs=TRUE, order_by="qstart") + theme_bw()
minimap2_long_ali_5000 <- dotplot(long_ali_5000, label_seqs=TRUE, order_by="qstart") + theme_bw()

# look at alignments with and without filtering
minimap2_aln 
minimap2_aln + minimap2_long_ali_5000

# see lots of qname chr which have just a few seqs
# lets focus on Z chr
W_vs_Zs <- all_ali[all_ali$tname %in% c('SUPER_Z1', 'SUPER_Z2'),]
minimap2_aln <- dotplot(W_vs_Zs, label_seqs=TRUE, order_by="qstart") + theme_bw()
minimap2_aln # raw dotplot without filtering


W_vs_Zs$qname

# de = Gap-compressed per-base sequence divergence
# multiply de by alignment length for each alignment
all_ali$de_scaled_alen <- all_ali$de * all_ali$alen
W1_vs_Z1 <- all_ali[(all_ali$qname == 'SUPER_W1') & (all_ali$tname == 'SUPER_Z1'),]
W1_vs_Z2 <- all_ali[(all_ali$qname == 'SUPER_W1') & (all_ali$tname == 'SUPER_Z2'),]
W2_vs_Z2 <- all_ali[(all_ali$qname == 'SUPER_W2') & (all_ali$tname == 'SUPER_Z2'),]

W1_Z1_nuc_divergence <- sum(W1_vs_Z1$de_scaled_alen) / sum(W1_vs_Z1$alen)
W1_Z2_nuc_divergence <- sum(W1_vs_Z2$de_scaled_alen) / sum(W1_vs_Z2$alen)
W2_Z2_nuc_divergence <- sum(W2_vs_Z2$de_scaled_alen) / sum(W2_vs_Z2$alen)

print(paste('Divergence for W1 Z1:', W1_Z1_nuc_divergence))
print(paste('Divergence for W1 Z2:', W1_Z2_nuc_divergence))
print(paste('Divergence for W2 Z2:', W2_Z2_nuc_divergence))

# percent of W1 covered by Z1
# The alen tag gives you the actual number of aligned bases (excluding gaps), 
# while the qstop - qstart difference provides the total alignment span on the query sequence, 
# which includes any gaps or unaligned regions.
# lets use alen to get percentage covered by an alignment
per_W1_aligned_by_Z1 <- (sum(W1_vs_Z1$alen)/W1_vs_Z1$qlen[1])*100
per_W1_aligned_by_Z2 <- (sum(W1_vs_Z2$alen)/W1_vs_Z2$qlen[1])*100
per_W2_aligned_by_Z2 <- (sum(W2_vs_Z2$alen)/W2_vs_Z2$qlen[1])*100

print(paste('Aln span for W1 Z1:', per_W1_aligned_by_Z1))
print(paste('Aln span for W1 Z2:', per_W1_aligned_by_Z2))
print(paste('Aln span for W2 Z2:', per_W2_aligned_by_Z2))


#W1_vs_Z1$q_aln <- W1_vs_Z1$qend - W1_vs_Z1$qstart
#(sum(W1_vs_Z1$q_aln)/W1_vs_Z1$qlen[1])*100

# Summary:
# 2.4% for W1 vs Z1
# 3.9% for W1 vs Z2
# 1.2% for W2 vs Z2


# filter alignments obtained from doing whole genome vs Ws to just get the Zs
long_ali_5000_filt_Zs <- long_ali_5000[long_ali_5000$tname %in% c('SUPER_Z1', 'SUPER_Z2'),]
minimap2_all_long_ali_5000_filt_Zs <- dotplot(long_ali_5000_filt_Zs, label_seqs=TRUE, order_by="qstart") + theme_bw()

minimap2_all_long_ali_5000 <- dotplot(all_long_ali_5000, label_seqs=TRUE, order_by="qstart") + theme_bw()
minimap2_all_long_ali_5000

plot_coverage(prim_alignment, fill='qname') + scale_fill_brewer(palette="Set1")

# plot on W chr what the coverage looks like
#all_long_ali_filt <- all_long_ali_5000[all_long_ali_5000$tname %in% c('SUPER_Z1', 'SUPER_Z2', 'SUPER_203', 'SUPER_35'),]
all_long_ali_filt <- all_long_ali_5000[all_long_ali_5000$tname %in% c('SUPER_Z1', 'SUPER_Z2'),]

all_flipped_long_ali_filt <- all_flipped_long_ali_5000[all_flipped_long_ali_5000$qname %in% c('SUPER_Z1', 'SUPER_Z2'),]
cov_plot_flipped <- plot_coverage(all_flipped_long_ali_filt, fill='qname')  + scale_fill_manual(values=c('#ff006e', '#8338ec'))
cov_plot_flipped
cov_plot <- plot_coverage(all_long_ali_filt, fill='qname')  + scale_fill_manual(values=c('#ffbe0b', '#fb5607'))
cov_plot


# Save plots
ggsave(plot=minimap2_all_long_ali_5000_filt_Zs, filename = "minimap2_SUPERWs_vs_SUPERZs_from_whole_genome_long_alignments_5000.pdf")
ggsave(plot=minimap2_aln, filename = "minimap2_SUPERWs_vs_SUPERZs.pdf")
ggsave(plot=minimap2_long_ali_10000, filename = "minimap2_SUPERWs_vs_SUPERZs_long_alignments_10000.pdf")
ggsave(plot=minimap2_long_ali_5000, filename = "minimap2_SUPERWs_vs_SUPERZs_long_alignments_5000.pdf")
ggsave(plot=cov_plot_flipped, filename = "coverage_plot_minimap2_SUPERWs_vs_all_long_alignments_5000_filt.pdf", height=4, width=6)
ggsave(plot=cov_plot, filename = "coverage_plot_minimap2_SUPERZs_from_Ws_vs_all_long_alignments_5000_filt.pdf", height=4, width=6)


# get coverage of Ws per chr
# Calculate the total alignment length per query chromosome and target chromosome

# get lengths of matches
all_long_ali_5000$t_length <- abs(all_long_ali_5000$tend - all_long_ali_5000$tstart)
all_long_ali_5000_W1 <- all_long_ali_5000[all_long_ali_5000$tname == 'SUPER_W1',]
all_long_ali_5000_W2 <- all_long_ali_5000[all_long_ali_5000$tname == 'SUPER_W2',]

alignment_lengths_W1 <- all_long_ali_5000_W1 %>%
  group_by(qname, tname) %>%
  summarise(sum_alns = sum(t_length))


alignment_lengths_W2 <- all_long_ali_5000_W2 %>%
  group_by(qname, tname) %>%
  summarise(sum_alns = sum(t_length))

alignment_lengths_W1$SUPER_W1_length <- 65809200
alignment_lengths_W2$SUPER_W2_length <- 4527220

alignment_lengths_W1$prop_covered <- (alignment_lengths_W1$sum_alns / alignment_lengths_W1$SUPER_W1_length)*100
alignment_lengths_W2$prop_covered <- (alignment_lengths_W2$sum_alns / alignment_lengths_W2$SUPER_W2_length)*100



