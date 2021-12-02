library(ggplot2)
library(forcats)

'%!in%' <- function(x,y)!('%in%'(x,y))

df <- read.csv('seq_depth.log', sep='\t', header=F)
meta <- read.csv('metadata.filtered.csv', header=T, row.names = 1)

disease.old <- c("Healthy","ACVD","advanced-adenoma","Crohns-disease","CRC",
                 "Obesity","Overweight","Rheumatoid-Arthritis","T2D", "Ulcerative-colitis")
disease.new <- c("Healthy","ACVD","AA","CD","CRC",
                 "OB","OW","RA","T2D", "UC")

test <- merge(df,meta,by.x='V1',by.y='SRA')

test <- test[test$disease %!in% c('Symptomatic-atherosclerosis','IGT','Underweight'),]

g <- ggplot(test, aes(x=disease, y=V6, group=disease, fill=disease))+
  geom_boxplot()+
  theme_bw()+
  ylab('Filtered Read Depth (RPM)')+
  xlab('Disesase')+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

# reorder
args <- c(list(g$data$disease), setNames(disease.old,disease.new))
g$data$disease <- factor(do.call(fct_recode,args), levels=disease.new)

ggsave(filename = 'read_depth.png',
       width = 12,
       height=6,
       units = 'in')
