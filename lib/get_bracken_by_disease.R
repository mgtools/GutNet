library(phyloseq)
library(microbiome)
library(vegan)
library(stringr)
suppressPackageStartupMessages(library("argparse"))

#setwd('~/Desktop/Microbiome_Network/EDA/')

####################################
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-i", "--input", help="Input Bracken output BIOM format")
parser$add_argument("-m", "--metadata", help="SRA metadata of healthy/disease phenotype")
parser$add_argument("-o", "--outdirr", help="Output directory to save split and filtered matrices")
parser$add_argument("--ignore", help="list of SRA's to ignore")
parser$add_argument("--filterdepth", default=1000000, type="double",
                    help="Sequencing depth threshold to filter samples at [default %(default)s]")
parser$add_argument("--prevelance filter", default=0.5, type="double",
                    help="Keep species with prevelance greater than, discard less than [default %(default)s]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# declare not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# import biom
ps <- import_biom(args$input, parseFunction=parse_taxonomy_default)

# tax name
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species")

# sra samples to remove
filter <- as.list(read.csv(args$ignore, header=F)[,1])

# import metadata
metadata <- read.csv(args$metadata,row.names=1)
row.names(metadata) <- metadata$SRA
metadata$Healthy <- 'Healthy'
metadata$Healthy[metadata$disease != 'Healthy'] <- 'Diseased'

# rename samples
sample_names(ps) <- str_remove(string = sample_names(ps),pattern = ".s")

# import metadata to phylo
sample_data(ps) <- metadata

# filter on read depth 1Mil
toRemove <- as.list(sample_names(ps)[colSums(otu_table(ps)) < args$filterdepth])

# merge filter and toRemove 
toRemove <- as.character(union(toRemove,filter))

# remove filtered samples
ps <- prune_samples(!(sample_names(ps) %in% toRemove), ps)

# remove SA samples
SA <- (sample_data(ps)[sample_data(ps)[,2] == 'Symptomatic-atherosclerosis'])$SRA
ps <- prune_samples(!(sample_names(ps) %in% SA), ps)

# get relative abundance 
#ps.relabu <- transform_sample_counts(ps, function(x) x / sum(x) )

# remove low prevelance taxa
filter_prevelance <- function(x){
  x.f <- prune_taxa((rowSums(otu_table(x) > 0) / ncol(otu_table(x))) >= 0.5,x) # get proportion of samples with > 0; remove less than 50%
  return(x.f)
}

# subset samples by disease
AA <- subset_samples(ps, disease == 'advanced-adenoma')
ACVD <- subset_samples(ps, disease == 'ACVD')
CRC <- subset_samples(ps, disease == 'CRC')
CD <- subset_samples(ps, disease == 'Crohns-disease')
WT <- subset_samples(ps, disease == 'Healthy')
IGT <- subset_samples(ps, disease == 'IGT')
OB <- subset_samples(ps, disease == 'Obesity')
OW <- subset_samples(ps, disease == 'Overweight')
RA <- subset_samples(ps, disease == 'Rheumatoid-Arthritis')
T2D <- subset_samples(ps, disease == 'T2D')
UC <- subset_samples(ps, disease == 'Ulcerative-colitis')
UW <- subset_samples(ps, disease == 'Underweight')

## filter by prevelance
AA <- filter_prevelance(AA)
ACVD <- filter_prevelance(ACVD)
CRC <- filter_prevelance(CRC)
CD <- filter_prevelance(CD)
WT <- filter_prevelance(WT)
IGT <- filter_prevelance(IGT)
OB <- filter_prevelance(OB)
OW <- filter_prevelance(OW)
RA <- filter_prevelance(RA)
T2D <- filter_prevelance(T2D)
UC <- filter_prevelance(UC)
UW <- filter_prevelance(UW)

# agglomerate to genus level
AA.g <- tax_glom(AA, taxrank='Genus', NArm=TRUE)
ACVD.g <- tax_glom(ACVD, taxrank='Genus', NArm=TRUE)
CRC.g <- tax_glom(CRC, taxrank='Genus', NArm=TRUE)
CD.g <- tax_glom(CD, taxrank='Genus', NArm=TRUE)
WT.g <- tax_glom(WT, taxrank='Genus', NArm=TRUE)
IGT.g <- tax_glom(IGT, taxrank='Genus', NArm=TRUE)
OB.g <- tax_glom(OB, taxrank='Genus', NArm=TRUE)
OW.g <- tax_glom(OW, taxrank='Genus', NArm=TRUE)
RA.g <- tax_glom(RA, taxrank='Genus', NArm=TRUE)
T2D.g <- tax_glom(T2D, taxrank='Genus', NArm=TRUE)
UC.g <- tax_glom(UC, taxrank='Genus', NArm=TRUE)
UW.g <- tax_glom(UW, taxrank='Genus', NArm=TRUE)

make_spp_table <- function(disease_ps){
  # makt otu table
  mtx_ps <- as.data.frame(otu_table(disease_ps))
  
  # make tax table
  taxmat <- as.data.frame(tax_table(disease_ps))
  
  # get spp name
  taxmat$Genus<-stringr::str_split_fixed(taxmat$Genus,'__',2)[,2]
  taxmat$Species<-stringr::str_split_fixed(taxmat$Species,'__',2)[,2]
  taxmat$name <- paste(taxmat$Genus, taxmat$Species)
  
  # rename rows
  row.names(mtx_ps) <- taxmat[row.names(mtx_ps),'name']
  
  return(mtx_ps)
}

make_genus_table <- function(disease_ps){
  # makt otu table
  mtx_ps <- as.data.frame(otu_table(disease_ps))
  
  # make tax table
  taxmat <- as.data.frame(tax_table(disease_ps))
  
  # rename genus
  taxmat$Genus<-stringr::str_split_fixed(taxmat$Genus,'__',2)[,2]
  
  # make genus unique
  taxmat$Genus <- make.unique(taxmat$Genus)
  
  # rename rows
  row.names(mtx_ps) <- taxmat[row.names(mtx_ps),'Genus']
  
  return(mtx_ps)
}

# make species otu table
AA_mtx <- make_spp_table(AA)
ACVD_mtx <- make_spp_table(ACVD)
CRC_mtx <- make_spp_table(CRC)
CD_mtx <- make_spp_table(CD)
WT_mtx <- make_spp_table(WT)
IGT_mtx <- make_spp_table(IGT)
OB_mtx <- make_spp_table(OB)
OW_mtx <- make_spp_table(OW)
RA_mtx <- make_spp_table(RA)
T2D_mtx <- make_spp_table(T2D)
UC_mtx <- make_spp_table(UC)
UW_mtx <- make_spp_table(UW)


# make species otu table
AA_mtx.g <- make_genus_table(AA.g)
ACVD_mtx.g <- make_genus_table(ACVD.g)
CRC_mtx.g <- make_genus_table(CRC.g)
CD_mtx.g <- make_genus_table(CD.g)
WT_mtx.g <- make_genus_table(WT.g)
IGT_mtx.g <- make_genus_table(IGT.g)
OB_mtx.g <- make_genus_table(OB.g)
OW_mtx.g <- make_genus_table(OW.g)
RA_mtx.g <- make_genus_table(RA.g)
T2D_mtx.g <- make_genus_table(T2D.g)
UC_mtx.g <- make_genus_table(UC.g)
UW_mtx.g <- make_genus_table(UW.g)

# save Rdata object
save.image(paste0(args$outdirr,'/Bracken_phyloseq_disease.filtered.RData'))

# save species csv 
write.csv(x = AA_mtx, file = paste0(args$outdirr, '/AA_bracken_filtered.s.csv'))
write.csv(x = ACVD_mtx, file = paste0(args$outdirr, '/ACVD_bracken_filtered.s.csv'))
write.csv(x = CRC_mtx, file = paste0(args$outdirr, '/CRC_bracken_filtered.s.csv'))
write.csv(x = CD_mtx, file = paste0(args$outdirr, '/CD_bracken_filtered.s.csv'))
write.csv(x = WT_mtx, file = paste0(args$outdirr, '/WT_bracken_filtered.s.csv'))
write.csv(x = IGT_mtx, file = paste0(args$outdirr, '/IGT_bracken_filtered.s.csv'))
write.csv(x = OB_mtx, file = paste0(args$outdirr, '/OB_bracken_filtered.s.csv'))
write.csv(x = OW_mtx, file = paste0(args$outdirr, '/OW_bracken_filtered.s.csv'))
write.csv(x = RA_mtx, file = paste0(args$outdirr, '/RA_bracken_filtered.s.csv'))
write.csv(x = T2D_mtx, file = paste0(args$outdirr, '/T2D_bracken_filtered.s.csv'))
write.csv(x = UC_mtx, file = paste0(args$outdirr, '/UC_bracken_filtered.s.csv'))
write.csv(x = UW_mtx, file = paste0(args$outdirr, '/UW_bracken_filtered.s.csv'))


# save species csv 
write.csv(x = AA_mtx.g, file = paste0(args$outdirr, '/AA_bracken_filtered.g.csv'))
write.csv(x = ACVD_mtx.g, file = paste0(args$outdirr, '/ACVD_bracken_filtered.g.csv'))
write.csv(x = CRC_mtx.g, file = paste0(args$outdirr, '/CRC_bracken_filtered.g.csv'))
write.csv(x = CD_mtx.g, file = paste0(args$outdirr, '/CD_bracken_filtered.g.csv'))
write.csv(x = WT_mtx.g, file = paste0(args$outdirr, '/WT_bracken_filtered.g.csv'))
write.csv(x = IGT_mtx.g, file = paste0(args$outdirr, '/IGT_bracken_filtered.g.csv'))
write.csv(x = OB_mtx.g, file = paste0(args$outdirr, '/OB_bracken_filtered.g.csv'))
write.csv(x = OW_mtx.g, file = paste0(args$outdirr, '/OW_bracken_filtered.g.csv'))
write.csv(x = RA_mtx.g, file = paste0(args$outdirr, '/RA_bracken_filtered.g.csv'))
write.csv(x = T2D_mtx.g, file = paste0(args$outdirr, '/T2D_bracken_filtered.g.csv'))
write.csv(x = UC_mtx.g, file = paste0(args$outdirr, '/UC_bracken_filtered.g.csv'))
write.csv(x = UW_mtx.g, file = paste0(args$outdirr, '/UW_bracken_filtered.g.csv'))
