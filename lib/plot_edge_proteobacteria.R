library(ggplot2)
library(reshape2)
library(igraph)
library(dplyr)
library(ggpubr)
library("viridis") 
library(reshape2)
library(ggrepel)

load('Networks_w_attributes.Rdata')

get_edge <- function(g){
  # get edge information
  edge_list <- get.data.frame(g, what="edges")
  if (ncol(edge_list) == 6){
    colnames(edge_list) <- c('head','tail','weight', 'weight2','sign','sign1')
  } else {
    colnames(edge_list) <- c('head','tail','weight','sign')
  }
  
  # get node attributes
  node_attributes <- data.frame(label = V(g)$label,
                                row.names   = V(g)$name,
                                identity = V(g)$identity,
                                degree      = V(g)$degree,
                                #closeness   = V(g)$closeness,
                                betweenness = V(g)$betweenness,
                                eigenvector = V(g)$eig,
                                Kingdom = V(g)$Kingdom,
                                Phylum = V(g)$Phylum,
                                Class = V(g)$Class,
                                Order = V(g)$Order,
                                Family = V(g)$Family,
                                Genus = V(g)$Genus)
  
  # populate taxa information
  edge_list$head_identity <- node_attributes[edge_list$head, 'identity']
  edge_list$tail_identity <- node_attributes[edge_list$tail, 'identity']
  edge_list$head_Kingdom <- node_attributes[edge_list$head, 'Kingdom']
  edge_list$tail_Kingdom <- node_attributes[edge_list$tail, 'Kingdom']
  edge_list$head_Phylum <- node_attributes[edge_list$head, 'Phylum']
  edge_list$tail_Phylum <- node_attributes[edge_list$tail, 'Phylum']
  edge_list$head_Class <- node_attributes[edge_list$head, 'Class']
  edge_list$tail_Class <- node_attributes[edge_list$tail, 'Class']
  edge_list$head_Order <- node_attributes[edge_list$head, 'Order']
  edge_list$tail_Order <- node_attributes[edge_list$tail, 'Order']
  edge_list$head_Family <- node_attributes[edge_list$head, 'Family']
  edge_list$tail_Family <- node_attributes[edge_list$tail, 'Family']
  edge_list$head_Genus <- node_attributes[edge_list$head, 'Genus']
  edge_list$tail_Genus <- node_attributes[edge_list$tail, 'Genus']
  
  count <- as.numeric(summary(edge_list$head_Phylum == 'Proteobacteria' | edge_list$tail_Phylum == 'Proteobacteria')[3])
  total <- nrow(edge_list)
  print(count/total)
  return(count/total)
}

WT_fw_s_proteobacteria <- get_edge(WT_fw.s)
WT_se_s_proteobacteria <- get_edge(WT_se.s)
AA_fw_s_proteobacteria <- get_edge(AA_fw.s)
AA_se_s_proteobacteria <- get_edge(AA_se.s)
ACVD_fw_s_proteobacteria <- get_edge(ACVD_fw.s)
ACVD_se_s_proteobacteria <- get_edge(ACVD_se.s)
CD_fw_s_proteobacteria <- get_edge(CD_fw.s)
CD_se_s_proteobacteria <- get_edge(CD_se.s)
CRC_fw_s_proteobacteria <- get_edge(CRC_fw.s)
CRC_se_s_proteobacteria <- get_edge(CRC_se.s)
OB_fw_s_proteobacteria <- get_edge(OB_fw.s)
OB_se_s_proteobacteria <- get_edge(OB_se.s)
OW_fw_s_proteobacteria <- get_edge(OW_fw.s)
OW_se_s_proteobacteria <- get_edge(OW_se.s)
RA_fw_s_proteobacteria <- get_edge(RA_fw.s)
RA_se_s_proteobacteria <- get_edge(RA_se.s)
T2D_fw_s_proteobacteria <- get_edge(T2D_fw.s)
T2D_se_s_proteobacteria <- get_edge(T2D_se.s)
UC_fw_s_proteobacteria <- get_edge(UC_fw.s)
UC_se_s_proteobacteria <- get_edge(UC_se.s)

listoflists <- list(WT_se_s_proteobacteria,
                    AA_se_s_proteobacteria,
                    ACVD_se_s_proteobacteria,
                    CD_se_s_proteobacteria,
                    CRC_se_s_proteobacteria,
                    OB_se_s_proteobacteria,
                    OW_se_s_proteobacteria,
                    RA_se_s_proteobacteria,
                    T2D_se_s_proteobacteria,
                    UC_se_s_proteobacteria
                    )
df <- as.data.frame(t(as.data.frame(listoflists)))
rownames(df)<-NULL
row.names(df) <- c('WT','AA','ACVD','CD','CRC','OB','OW','RA','T2D','UC')
colnames(df) <- c('SPIEC-Easi')
df$disease <- row.names(df)
df.melt <- melt(df)

df.melt$label <- NA
df.melt[df.melt$disease=='WT','label'] <- 'Healthy'

g <- ggplot(df.melt, aes(x=variable,y=value))+
  #geom_violin()+
  geom_boxplot(width=0.3)+
  #geom_jitter(width=0.1)+
  geom_label_repel(aes(label=label), direction='x',nudge_x = 0.3)+
  theme_bw()+
  scale_y_continuous(limits = c(0.25,0.75), labels = scales::percent)+
  ylab('Proportion of Edges w/Proteobacteria (%)')+
  #xlab('Association Network Method')+
  xlab('')+
  ggtitle('Proportion of Proteobacteria Edges', 'Healthy and Disease Networks')

g

ggsave('Edge_protobacteria_boxplot.png',g,height=6,width=5,units='in')
