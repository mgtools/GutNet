library(ggplot2)
library(dplyr)
library(ggpubr)
library(igraph)
library(patchwork)
library(forcats)
library(reshape2)

load('Bracken_phyloseq_disease.filtered.RData')
rm(list=setdiff(ls(), "ps"))


WT <- filter_taxa(subset_samples(ps, disease=="Healthy"), function(y) sum(y) > 0, TRUE)
AA <- filter_taxa(subset_samples(ps, disease=="advanced-adenoma"), function(y) sum(y) > 0, TRUE)
ACVD <- filter_taxa(subset_samples(ps, disease=="ACVD"), function(y) sum(y) > 0, TRUE)
CD <- filter_taxa(subset_samples(ps, disease=="Crohns-disease"), function(y) sum(y) > 0, TRUE)
CRC <- filter_taxa(subset_samples(ps, disease=="CRC"), function(y) sum(y) > 0, TRUE)
OB <- filter_taxa(subset_samples(ps, disease=="Obesity"), function(y) sum(y) > 0, TRUE)
OW <- filter_taxa(subset_samples(ps, disease=="Overweight"), function(y) sum(y) > 0, TRUE)
T2D <- filter_taxa(subset_samples(ps, disease=="T2D"), function(y) sum(y) > 0, TRUE)
RA <- filter_taxa(subset_samples(ps, disease=="Rheumatoid-Arthritis"), function(y) sum(y) > 0, TRUE)
UC <- filter_taxa(subset_samples(ps, disease=="Ulcerative-colitis"), function(y) sum(y) > 0, TRUE)

WT <- transform_sample_counts(WT, function(x) x / sum(x) )
AA <- transform_sample_counts(AA, function(x) x / sum(x) )
ACVD <- transform_sample_counts(ACVD, function(x) x / sum(x) )
CD <- transform_sample_counts(CD, function(x) x / sum(x) )
CRC <- transform_sample_counts(CRC, function(x) x / sum(x) )
OB <- transform_sample_counts(OB, function(x) x / sum(x) )
OW <- transform_sample_counts(OW, function(x) x / sum(x) )
T2D <- transform_sample_counts(T2D, function(x) x / sum(x) )
RA <- transform_sample_counts(RA, function(x) x / sum(x) )
UC <- transform_sample_counts(UC, function(x) x / sum(x) )

diseases = c('Healthy',"AA","ACVD","CD","CRC","OB", 'OW',"RA","T2D","UC")


df = data.frame(matrix(ncol=21,nrow=10, dimnames=list(diseases,as.character(seq(0,1,0.05)))), check.names=F)

# build function
get_ntax <- function(ps,x,df){

  for (i in seq(0,1,0.05)){
    if (i == 1){
      df[x,'0'] = 0 
      next
    }
    taxa_count = (ntaxa(filter_taxa(ps, function(x) sum(x > 0) > (i*length(x)), TRUE)))
    n = as.character(round(1-i, digits=2))
    #i = as.character(i)
    print(c(x,n,taxa_count))
    df[x,n] = taxa_count/ntaxa(ps)
    #print(i)
  }
  
  return(df)
}

# get dataframe
df = get_ntax(WT,'Healthy',df)
df = get_ntax(AA,'AA',df)
df = get_ntax(ACVD,'ACVD',df)
df = get_ntax(CD,'CD',df)
df = get_ntax(CRC,'CRC',df)
df = get_ntax(OB,'OB',df)
df = get_ntax(OW,'OW',df)
df = get_ntax(RA,'RA',df)
df = get_ntax(T2D,'T2D',df)
df = get_ntax(UC,'UC',df)

df <- within(df, rm('0'))

# melt
df.melt <- melt(as.matrix(df))

df.melt$Var1 <- factor(df.melt$Var1)

# plot
g <- ggplot(df.melt, aes(x=Var2,y=value,color=Var1,group=Var1))+
  scale_shape_manual(values=1:10) +
  geom_point(aes(shape=Var1))+
  geom_line()+
  scale_x_reverse(breaks = seq(0.05,1.00,0.05),labels = scales::percent)+
  #scale_x_continuous(breaks = seq(0,1,0.05),labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  ylab('Species Present (%)')+
  xlab('Prevelance (% of samples with species present)')

g

g$labels$colour <- "Diseases"
g$labels$shape <- "Diseases"

ggsave(filename = 'spp_prevelance_percent.png',
       plot = g, width = 12, height = 6, unit='in')





