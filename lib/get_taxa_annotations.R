suppressPackageStartupMessages(library(hash))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(reshape2))

# load df
ps <- import_biom('sung_bracken.biom', parseFunction=parse_taxonomy_default)

# tax name
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order",
                             "Family", "Genus", "Species")

# get tax table 
taxmat <- as.data.frame(tax_table(ps))

# get spp name
taxmat$name <- paste(stringr::str_split_fixed(taxmat$Genus,'__',2)[,2], stringr::str_split_fixed(taxmat$Species,'__',2)[,2])

# get id
taxmat$id <- row.names(taxmat)

# reorder
taxmat <- taxmat[,c(9:8,1:7)]

write.csv(taxmat,'0_taxa_annotations/phyloseq_id_to_annotation.csv')


# get genus
ps.g <- tax_glom(ps, taxrank='Genus', NArm=TRUE)

# tax table
taxmat.g <- as.data.frame(tax_table(ps.g))

# rename genus
taxmat.g$name <-stringr::str_split_fixed(taxmat.g$Genus,'__',2)[,2]

# make genus unique
taxmat.g$name <- make.unique(taxmat.g$name)

# get id
taxmat.g$id <- row.names(taxmat.g)

# reorder
taxmat.g <- taxmat.g[,c(9:8,1:7)]

write.csv(taxmat.g,'0_taxa_annotations/phyloseq_id_to_annotation.g.csv')
