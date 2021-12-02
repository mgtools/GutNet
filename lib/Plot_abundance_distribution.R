library(ggplot2)
library(reshape2)
library(igraph)
library(dplyr)
library(ggpubr)
library(phyloseq)
library(microbiome)
library(forcats)

'%!in%' <- function(x,y)!('%in%'(x,y))
load('Bracken_phyloseq_disease.filtered.RData')

# agglomorate at Phylum
ps.p <- tax_glom(ps,taxrank='Phylum')

# remove IGT, SA, UW
ps.p.x <- prune_samples(sample_data(ps.p)$disease %!in% c('Underweight','IGT','SA'), ps.p)

# rel abu
ps.p.e <- transform_sample_counts(ps.p.x, function(OTU) OTU/sum(OTU) )

# parse phylum information
abu <- (as.data.frame(otu_table(ps.p.e)))
tax_tab <- as.data.frame(tax_table(ps.p.e))
tax_tab$Phylum <- as.data.frame(stringr::str_split_fixed(tax_tab$Phylum,'__',2))[2]

# rename rows
row.names(abu) <- tax_tab[row.names(abu),2]$V2

# reorder by sum abundances
abu <- abu[rev(order(rowSums(abu),row.names(abu))),]

abu_1 <- abu[1:5,]
abu_2 <- t(colSums(abu[6:nrow(abu),]))
rownames(abu_2) <- c('Other')

abu_new <- rbind(abu_1,abu_2)

# melt
abu.melt <- melt(as.matrix(abu_new))

# get metadata
sd <- as.matrix(sample_data(ps.p.e))

# fill disease
abu.melt$disease <- sd[abu.melt$Var2,'disease']

# rename
colnames(abu.melt) <- c('Phylum','SRA','count','Phenotype')

cust_colors <- rev(c('#B8BDB5','#005F73','#0A9396','#94D2BD','#E9D8A6',
                   '#EE9B00'))

# rename 
disease.old <- c("Healthy","ACVD","advanced-adenoma","Crohns-disease","CRC",
                 "Obesity","Overweight","Rheumatoid-Arthritis","T2D", "Ulcerative-colitis")
disease.new <- c("Healthy","ACVD","AA","CD","CRC",
                 "OB","OW","RA","T2D", "UC")

args <- c(list(abu.melt$Phenotype), setNames(disease.old,disease.new))
abu.melt$Phenotype <- factor(do.call(fct_recode,args), levels=disease.new)

# plot distribution of phylum
g <- ggplot(abu.melt, aes(x=Phylum, y=count,fill=factor(Phylum)))+
  geom_boxplot()+
  facet_wrap(~Phenotype,ncol = 5,nrow = 2)+
  scale_fill_manual(values=cust_colors,aesthetics = c("colour", "fill"))+ 
  ylab('Relative Abundance (%)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('~/Desktop/Microbiome_Network/EDA/Phylumn_abundance_distribution_boxplot.png',
       plot = g, 
       height=8, 
       width=16, 
       units='in')

# get means
abu.mean <- abu.melt %>%
  group_by(Phylum,Phenotype) %>%
  dplyr::summarize(Mean = mean(count, na.rm=TRUE),
                   Mean_p = percent(Mean, accuracy=0.1),
                   n = n(),
                   sd = sd(count),
                   std = sd/sqrt(n),
                   std_p = percent(std, accuracy=0.1),
                   lab = paste0(Mean_p,'\n','(','\u00B1',std_p,')'),
                   grouping = 'By Diseases')

# get all info
abu.all <- abu.melt %>%
  group_by(Phylum) %>%
  dplyr::summarize(Phenotype = 'All',
                   Mean = mean(count, na.rm=TRUE),
                   Mean_p = percent(Mean, accuracy=0.1),
                   n = n(),
                   sd = sd(count),
                   std = sd/sqrt(n),
                   std_p = percent(std, accuracy=0.1),
                   lab = paste0(Mean_p,'\n','(','\u00B1',std_p,')'),
                   grouping = 'All')

# merge
abu.all.merge <- rbind(abu.mean, abu.all)

#abu.mean$order <- rev(order(abu.mean$Mean))

# remove label for lower abu phylums 
abu.all.merge$lab <- ifelse((abu.all.merge$Phylum %in% c('Bacteroidetes','Firmicutes','Proteobacteria') | (abu.all.merge$Phylum == 'Actinobacteria' & abu.all.merge$Mean >0.05)), abu.all.merge$lab, '')


# plot
h <- ggplot(abu.all.merge, aes(x=Phenotype, y=Mean,fill=factor(Phylum)))+
  geom_bar(stat='identity', position="fill", colour="black")+
  scale_fill_manual(values=cust_colors,aesthetics = c("colour", "fill"))+
  geom_text(aes(label=lab),position=position_stack(0.5))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  ylab('Mean Relative Abundance (%)')+
  ggtitle('Mean Relative Abundance - Phylum level')+
  theme_bw()+
  theme(axis.title=element_text(size=15))+
  facet_grid(~grouping, space="free",scales="free")

h$labels$fill <- 'Phylum'
h
  
ggsave('Phylumn_abundance_distribution_average.png',
       plot = h, 
       height=8, 
       width=10.5, 
       units='in')
