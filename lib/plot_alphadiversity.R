 qlibrary(ggplot2)
library(microbiome)
library(tsnemicrobiota)
library(phyloseq)
library(gginnards)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(forcats)

'%!in%' <- function(x,y)!('%in%'(x,y))

load('Bracken_phyloseq_disease.filtered.RData')

# remove disease not analyzed
ps <- subset_samples(ps, disease%!in%c('IGT','Underweight','SA'))
# remove zero rows 
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

disease.old <- c("Healthy","ACVD","advanced-adenoma","Crohns-disease","CRC",
             "Obesity","Overweight","Rheumatoid-Arthritis","T2D", "Ulcerative-colitis")
disease.new <- c("Healthy","ACVD","AA","CD","CRC",
             "OB","OW","RA","T2D", "UC")

# plot all diseases 
r.disease.obs <- plot_richness(ps, x='disease' , measures=c("Observed"))+ 
          geom_violin(aes(fill=disease), trim = T, alpha=0.4) 
r.disease.obs$layers <- r.disease.obs$layers[-1]
r.disease.obs <- r.disease.obs+ 
          geom_boxplot(width=0.1,outlier.shape = NA)+
          theme_bw()+
          #stat_compare_means(method = "anova", label.y = 8000)+      # Add global p-value 
          stat_compare_means(ref.group='Healthy', 
                             method="wilcox.test", 
                             hide.ns = F, 
                             label = "..p.signif..", 
                             label.y=6300)+ 
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0, 6450))

# reorder
args <- c(list(r.disease.obs$data$disease), setNames(disease.old,disease.new))
r.disease.obs$data$disease <- factor(do.call(fct_recode,args), levels=disease.new)

# delete points
r.disease.obs <- delete_layers(r.disease.obs,'GeomPoint')

r.disease.shannon <- plot_richness(ps, x='disease' , measures=c("Shannon"))+ 
  geom_violin(aes(fill=disease), trim = T, alpha=0.4) 
r.disease.shannon$layers <- r.disease.shannon$layers[-1]
r.disease.shannon <- r.disease.shannon + 
  geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_bw()+
  #stat_compare_means(method = "anova", label.y = 7.5)+      # Add global p-value 
  stat_compare_means(ref.group='Healthy', 
                     method="wilcox.test", 
                     hide.ns = F, 
                     label = "..p.signif..", 
                     p.adjust.method="BH", 
                     label.y=6.5)+ 
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0, 7))
# reorder
args <- c(list(r.disease.shannon$data$disease), setNames(disease.old,disease.new))
r.disease.shannon$data$disease <- factor(do.call(fct_recode,args), levels=disease.new)

# delete points
r.disease.shannon <- delete_layers(r.disease.shannon,'GeomPoint')



# plot healthy vs diseased - obs
r.healthy.obs <- plot_richness(ps, x='Healthy', measures=c("Observed"))+ 
  geom_violin(aes(fill=Healthy), trim = T, alpha=0.4)
r.healthy.obs$layers <- r.healthy.obs$layers[-1]
r.healthy.obs <- r.healthy.obs+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_bw()+
  stat_compare_means(ref.group='Healthy', 
                     method="wilcox.test", 
                     hide.ns = F, 
                     label = "..p.signif..", 
                     label.y=6300)+ 
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0, 6450))
r.healthy.obs$data$Healthy <- factor(r.healthy.obs$data$Healthy, levels = c('Healthy','Diseased'))
r.healthy.obs <- r.healthy.obs + scale_fill_manual(values=c('#F8766D','#B8BDB5'))

# plot healthy vs disease - shannon
r.healthy.shannon <- plot_richness(ps, x='Healthy', measures=c("Shannon"))+ 
  geom_violin(aes(fill=Healthy), trim = T, alpha=0.4)
r.healthy.shannon$layers <- r.healthy.shannon$layers[-1]
r.healthy.shannon <- r.healthy.shannon+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  theme_bw()+
  stat_compare_means(ref.group='Healthy', 
                     method="wilcox.test", 
                     hide.ns = F, 
                     label = "..p.signif..", 
                     label.y=6.5)+ 
  theme(legend.position = 'none')+
  scale_y_continuous(limits = c(0, 7))
r.healthy.shannon$data$Healthy <- factor(r.healthy.shannon$data$Healthy, levels = c('Healthy','Diseased'))
r.healthy.shannon <- r.healthy.shannon + scale_fill_manual(values=c('#F8766D','#B8BDB5'))




# merge plots
#g <- r.disease.obs + r.disease.shannon & theme(legend.position = 'bottom')
#r.disease <- g+ plot_layout(guides = "collect")

g <- r.healthy.obs + xlab('') + ylab('Species Richness') + r.disease.obs + ylab('') + xlab('') + plot_layout(widths=c(1,4)) + r.healthy.shannon + xlab('') + ylab('Shannon diversity') + r.disease.shannon + ylab('') + xlab('') + plot_layout(widths=c(1,4))

#g <- ggarrange(observed, shannon, labels = c('A','B'), nrow=2, ncol=1)

g <- g + plot_annotation(tag_levels = 'A')

# save plot
ggsave('Alpha_Diversity_plots.png',
       plot = g,
       width = 8,
       height = 6,
       units = 'in')

