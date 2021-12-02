library(ggplot2)
library(microbiome)
library(tsnemicrobiota)
library(phyloseq)
library(ggpubr)
library(grid)
library(gridExtra)
library(patchwork)
load('~/Desktop/Microbiome_Network/Filtered_abundance_mtx/Bracken_phyloseq_disease.filtered.RData')

disease.old <- c("Healthy","ACVD","advanced-adenoma","Crohns-disease","CRC",
                 "Obesity","Overweight","Rheumatoid-Arthritis","T2D", "Ulcerative-colitis")
disease.new <- c("Healthy","ACVD","AA","CD","CRC",
                 "OB","OW","RA","T2D", "UC")

# remove disease not analyzed
ps <- subset_samples(ps, disease%!in%c('IGT','Underweight','SA'))
# remove zero rows 
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Convert to compositional data
ps.rel <- microbiome::transform(ps, "compositional")

# Pick core taxa with with the given prevalence and detection limits
ps.core <- core(ps.rel, detection = .1/100, prevalence = 90/100)

# Use relative abundances for the core
ps.core <- microbiome::transform(ps.core, "compositional")

# NMDS w/ brays - color by disease
ord <- ordinate(ps.core, "MDS", "bray")
g <- plot_ordination(ps.core, ord, color = "disease") +
  geom_point(size = 2)+
  ggtitle('NMDS (Bray-Curtis)')+
  theme_bw()

# reorder
args <- c(list(g$data$disease), setNames(disease.old,disease.new))
g$data$disease <- factor(do.call(fct_recode,args), levels=disease.new)

# colors
cust_colors <- c('#B8BDB5','#005F73','#0A9396','#94D2BD','#E9D8A6',
                 '#EE9B00','#CA6702','#BB3E03','#AE2012','#9B2226')

# color
g <- g + scale_color_manual(values = cust_colors)

# get ellipse
g.ellipse <- g+stat_ellipse()

# NMDS w/ brays - color by healthy/disease
h <- plot_ordination(ps.core, ord, color = "Healthy") +
  geom_point(size = 2, shape=16)+ggtitle('NMDS (Bray-Curtis)')+
  theme_bw()
h.ellipse <- h+stat_ellipse()

# plot tsne
tsne_res <- tsne_phyloseq(ps.core, distance="bray", perplexity = 30)

j<- plot_tsne_phyloseq(ps.core, tsne_res,
                   color = 'disease', title='t-SNE (Bray-Curtis)') +
  geom_point(size=2, shape=16)+
  theme_bw()

# reorder
args <- c(list(j$data$disease), setNames(disease.old,disease.new))
j$data$disease <- factor(do.call(fct_recode,args), levels=disease.new)

# colors
cust_colors <- c('#B8BDB5','#005F73','#0A9396','#94D2BD','#E9D8A6',
                 '#EE9B00','#CA6702','#BB3E03','#AE2012','#9B2226')

# color
j <- j + scale_color_manual(values = cust_colors)

j.ellipse <- j+stat_ellipse()

k<- plot_tsne_phyloseq(ps.core, tsne_res,
                       color = 'Healthy', title='t-SNE (Bray-Curtis)') +
  geom_point(size=2, shape=16)+
  theme_bw()
k.ellipse <- k+stat_ellipse()

# merge plots
combine <- g + j & theme(legend.position = 'bottom')
m <- combine+ plot_layout(guides = "collect")

combine <- h + k & theme(legend.position = 'bottom')
n <- combine+ plot_layout(guides = "collect")

ordination_plots <- ggarrange(n,m,nrow=2,ncol=1)

ggsave(filename = 'ordination_filtered_data.png',
       plot = ordination_plots,
       width = 12,
       height=12,
       units = 'in')


# merge plots - ellipse
combine <- g.ellipse + j.ellipse & theme(legend.position = 'bottom')
m <- combine+ plot_layout(guides = "collect")

combine <- h.ellipse + k.ellipse & theme(legend.position = 'bottom')
n <- combine+ plot_layout(guides = "collect")

ordination_plots <- ggarrange(n,m,nrow=2,ncol=1)

ggsave(filename = 'ordination_filtered_data.ellipse.png',
       plot = ordination_plots,
       width = 12,
       height=12,
       units = 'in')
