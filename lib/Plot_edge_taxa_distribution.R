library(ggplot2)
library(reshape2)
library(igraph)
library(dplyr)
library(ggpubr)
library("viridis") 

load('Networks_w_attributes.Rdata')
'%!in%' <- function(x,y)!('%in%'(x,y))
# determine top interaction type - positive and negative edges 

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

get_edge_by_taxa <- function(g,disease,taxa){
  
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
  
  counts_pos = data.frame()
  counts_neg = data.frame()
  
  # count anno information
  for (i in 1:nrow(edge_list)){
    h.anno <- edge_list[i,paste('head',taxa,sep = '_')]
    t.anno <- edge_list[i,paste('tail',taxa,sep = '_')]
    sign <- edge_list[i,'sign']
    
    # remove empty anno
    if (h.anno == "" || t.anno == ""){
      next
    }
    
    # populate if missing positive 
    if (sign == 'positive'){
      if (h.anno %!in% row.names(counts_pos)){
        counts_pos[h.anno,h.anno] = 0
      }
      if (t.anno %!in% row.names(counts_pos)){
        counts_pos[t.anno,t.anno] = 0
      }
      counts_pos[is.na(counts_pos)] = 0
    } else { # populate missing negative
      if (h.anno %!in% row.names(counts_neg)){
        counts_neg[h.anno,h.anno] = 0
      }
      if (t.anno %!in% row.names(counts_neg)){
        counts_neg[t.anno,t.anno] = 0
      }
      counts_neg[is.na(counts_neg)] = 0
    }
    
    # add count positive
    if (sign == 'positive'){
      if (h.anno == t.anno){
        counts_pos[h.anno,t.anno] = counts_pos[h.anno,t.anno] + 1
        next
      }
      counts_pos[h.anno,t.anno] = counts_pos[h.anno,t.anno] + 1
      counts_pos[t.anno,h.anno] = counts_pos[t.anno,h.anno] + 1
    } else {
      if (h.anno == t.anno){
        counts_neg[h.anno,t.anno] = counts_neg[h.anno,t.anno] + 1
        next
      }
      counts_neg[h.anno,t.anno] = counts_neg[h.anno,t.anno] + 1
      counts_neg[t.anno,h.anno] = counts_neg[t.anno,h.anno] + 1
    }
  }
  
  # fill zeros
  counts_pos[is.na(counts_pos)] = 0
  counts_neg[is.na(counts_neg)] = 0
  
  # fill lower tri with NA
  counts_pos[lower.tri(counts_pos,diag=F)] <- NA
  counts_pos.melt <- reshape2::melt(as.matrix(counts_pos), na.rm=T)
  counts_neg[lower.tri(counts_neg,diag=F)] <- NA
  counts_neg.melt <- reshape2::melt(as.matrix(counts_neg), na.rm=T)
  
  # rename
  colnames(counts_pos.melt) <- c('Head','Tail','Count')
  colnames(counts_neg.melt) <- c('Head','Tail','Count')
  
  # reorder 
  for (i in 1:nrow(counts_pos.melt)){
    counts_pos.melt[i,1:2] <- counts_pos.melt[i,2:1][order(counts_pos.melt[i,1:2])]
  }
  for (i in 1:nrow(counts_neg.melt)){
    counts_neg.melt[i,1:2] <- counts_neg.melt[i,2:1][order(counts_neg.melt[i,1:2])]
  }
  
  # add edge
  counts_pos.melt$edge <- paste(counts_pos.melt$Head, counts_pos.melt$Tail, sep='-')
  counts_pos.melt$disease <- disease
  counts_neg.melt$edge <- paste(counts_neg.melt$Head, counts_neg.melt$Tail, sep='-')
  counts_neg.melt$disease <- disease
  
  # add sign
  counts_pos.melt$sign <- "Positive"
  counts_neg.melt$sign <- "Negative"
  
  # merge
  counts.melt <- rbind(counts_pos.melt,counts_neg.melt)
  
  return(counts.melt)
}

taxa = 'Phylum'

AA_se.s.edgecounts <- get_edge_by_taxa(AA_se.s,'AA', taxa)
ACVD_se.s.edgecounts <- get_edge_by_taxa(ACVD_se.s,'ACVD', taxa)
CD_se.s.edgecounts <- get_edge_by_taxa(CD_se.s,'CD', taxa)
CRC_se.s.edgecounts <- get_edge_by_taxa(CRC_se.s,'CRC', taxa)
OB_se.s.edgecounts <- get_edge_by_taxa(OB_se.s,'OB', taxa)
OW_se.s.edgecounts <- get_edge_by_taxa(OW_se.s,'OW', taxa)
RA_se.s.edgecounts <- get_edge_by_taxa(RA_se.s,'RA', taxa)
T2D_se.s.edgecounts <- get_edge_by_taxa(T2D_se.s,'T2D', taxa)
UC_se.s.edgecounts <- get_edge_by_taxa(UC_se.s,'UC', taxa)
WT_se.s.edgecounts <- get_edge_by_taxa(WT_se.s,'WT', taxa)

merge_edge_df <- function(method, taxa){
  df <- merge(merge(merge(merge(merge(merge(merge(merge(merge(
    eval(as.name(paste0("AA_",method,'.',taxa,'.edgecounts'))),
    eval(as.name(paste0("ACVD_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("CD_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("CRC_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("WT_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("OB_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("OW_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("RA_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("T2D_",method,'.',taxa,'.edgecounts'))), all = TRUE),
    eval(as.name(paste0("UC_",method,'.',taxa,'.edgecounts'))), all = TRUE)
  return (df)
}

plot_stacked_bar <- function(df, taxa, n=9){
  
  # find top n taxa to keep, regroup everything else as other
  keep_pos <- df[df$sign=='Positive',] %>% 
    arrange(desc(Count)) %>% 
    group_by(disease) %>% slice(1:n) %>%
    pull(edge) %>%
    unique
  
  keep_neg <- df[df$sign=='Negative',] %>% 
    arrange(desc(Count)) %>% 
    group_by(disease) %>% slice(1:n) %>%
    pull(edge) %>%
    unique
  
  # colors
  cust_colors <- c(c('#B8BDB5','#005F73','#0A9396','#94D2BD','#E9D8A6',
                     '#EE9B00','#CA6702','#BB3E03','#AE2012','#9B2226'),sample(color, (length(keep_pos)+length(keep_neg)-9)))
  #cust_colors <- colors[1:length(keep_pos)+length(keep_neg)]
  
  
  # rename else to Other
  df$plot_lab <- NA
  df[df$sign=='Positive',]$plot_lab <- ifelse(df[df$sign=='Positive',]$edge %in% keep_pos, df[df$sign=='Positive',]$edge, "Other")
  df[df$sign=='Negative',]$plot_lab <- ifelse(df[df$sign=='Negative',]$edge %in% keep_neg, df[df$sign=='Negative',]$edge, "Other")
  
  # get factors
  lab <- sort(unique(df$plot_lab))
  lab <- lab[lab%!in%c('Other')]
  lab <- c('Other',lab)
  
  # get factors
  lab <- sort(unique(df$plot_lab))
  lab <- lab[lab%!in%c('Other')]
  lab <- c('Other',lab)
  
  # set disease order
  df$disease <- factor(df$disease, levels=c('WT',"AA","ACVD","CD","CRC",
                                            "OB","OW","RA","T2D", "UC"))
  
  # aggregate top n
  df.agg <- aggregate(Count ~ plot_lab + disease + sign + method, data=df, FUN=sum)
  
  # calculate totals 
  totals <- df.agg %>% group_by(disease,sign) %>% summarize(total=sum(Count))
  totals$lab <- paste0("n=",totals$total)
  
  # plot
  g <-ggplot(df.agg, aes(x=disease, y=Count, fill=factor(plot_lab, levels=lab)))+ 
    geom_bar(stat='identity', position="fill", colour="black")+
    scale_fill_manual(values=cust_colors,aesthetics = c("colour", "fill"))+
    ylab('Proportion (%)')+
    geom_text(data = totals, aes(x=disease, label=lab,fill=NULL,y=1,vjust=-1),size=2)+
    facet_grid(~sign)
  
  # chance label
  g$labels$fill <- taxa
  
  return(g)
}

plot_stacked_bar_fw_se <- function(df, taxa, n=9){
  
  keep_pos_se <- df[df$sign=='Positive' & df$method=='SPIEC-EASI' ,] %>% 
    arrange(desc(Count)) %>% 
    group_by(disease) %>% slice(1:n) %>%
    pull(edge) %>%
    unique
  
  keep_neg_se <- df[df$sign=='Negative'& df$method=='SPIEC-EASI' ,] %>% 
    arrange(desc(Count)) %>% 
    group_by(disease) %>% slice(1:n) %>%
    pull(edge) %>%
    unique
  
  keep <- unique(c(keep_pos_se,keep_neg_se))
  
  # colors
  cust_colors <- c(c('#B8BDB5','#005F73','#0A9396','#94D2BD','#E9D8A6',
                     '#EE9B00','#CA6702','#BB3E03'), ####'#AE2012','#9B2226'),
                      sample(color, (length(keep)+2)))
  
  # rename else to Other
  df$plot_lab <- NA
  df[df$sign=='Positive' & df$method=='SPIEC-EASI',]$plot_lab <- ifelse(df[df$sign=='Positive' & df$method=='SPIEC-EASI',]$edge %in% keep_pos_se, df[df$sign=='Positive' & df$method=='SPIEC-EASI',]$edge, "Other")
  df[df$sign=='Negative' & df$method=='SPIEC-EASI',]$plot_lab <- ifelse(df[df$sign=='Negative' & df$method=='SPIEC-EASI',]$edge %in% keep_neg_se, df[df$sign=='Negative' & df$method=='SPIEC-EASI',]$edge, "Other")
  
  # get factors
  lab <- sort(unique(df$plot_lab))
  lab <- lab[lab%!in%c('Other')]
  lab <- c('Other',lab)
  
  # get factors
  lab <- sort(unique(df$plot_lab))
  lab <- lab[lab%!in%c('Other')]
  lab <- c('Other',lab)
  
  # set disease order
  df$disease <- factor(df$disease, levels=c('WT',"AA","ACVD","CD","CRC",
                                            "OB","OW","RA","T2D", "UC"))
  
  # aggregate top n
  df.agg <- aggregate(Count ~ plot_lab + disease + sign + method, data=df, FUN=sum)
  
  # calculate totals 
  totals <- df.agg %>% group_by(disease,sign,method) %>% summarize(total=sum(Count))
  totals$lab <- paste0("n=",totals$total)
  
  # plot
  g <-ggplot(df.agg, aes(x=disease, y=Count, fill=factor(plot_lab, levels=lab)))+ 
    geom_bar(stat='identity', position="fill", colour="black")+
    scale_fill_manual(values=cust_colors,aesthetics = c("colour", "fill"))+
    ylab('Proportion (%)')+
    geom_text(data = totals, aes(x=disease, label=lab,fill=NULL,y=1,vjust=-1),size=2)+
    facet_wrap(method~sign,nrow = 2,ncol=2)
  
  # chance label
  g$labels$fill <- paste(taxa, taxa, sep='-')
  
  return(g)
}

# merge disease dataframe
edgecounts.se.s <- merge_edge_df('se','s')

# add annotations for method
edgecounts.se.s$method <- 'SPIEC-EASI'
edgecounts.se.s$sign <- factor(edgecounts.se.s$sign, levels=c('Positive','Negative'))

# remove zero counts
edgecounts.se.s <- edgecounts.se.s[edgecounts.se.s$Count != 0,]

# plot
g.s <- plot_stacked_bar_fw_se(edgecounts.s, taxa=taxa, n = 10)+ 
          ggtitle('Edge Phylum level Distribution (species-level network) - edge agglomerated @ Phylumn')+ 
          theme(legend.position = "bottom")

# save
ggsave('Edge_taxa_distribution_Phylum.s.png',plot = g.s, height=12, width=16, units='in')



