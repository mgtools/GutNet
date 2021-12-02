library(SpiecEasi)
library(Matrix)

args <- commandArgs(TRUE)

infile <- args[1]
output_image <- args[2]
output_mb_adj <- args[3]

print(infile)

# load dataframe
df <- read.csv(infile, row.names=1,header=TRUE)

# set dataframe as numerical matrix + transform
df.mtx <- t(as.matrix(df))
cn=colnames(df.mtx)
#print(df.mtx)
# set seed
set.seed(1)

# run spiec

print('running mb')
se.mb <- spiec.easi(df.mtx, method='mb', lambda.min.ratio=9e-5, nlambda=200, pulsar.params=list(rep.num=30, ncores=48))

# get mb edge weights
se.mb.beta <- symBeta(getOptBeta(se.mb),mode='maxabs')

# summary
se.mb.beta.summary <- summary(se.mb.beta)

# annotate
se.mb.beta.summary$i_annotated <- cn[se.mb.beta.summary$i]
se.mb.beta.summary$j_annotated <- cn[se.mb.beta.summary$j]

save.image(file=output_image)

# output adj list 
write.csv(se.mb.beta.summary, output_mb_adj)

