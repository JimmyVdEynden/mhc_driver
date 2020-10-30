###########################################
# Download cosmic data
###########################################

# CGC
######

# v91, 17/04/2020
# Manual download from https://cancer.sanger.ac.uk/cosmic/download

# General data (no difference tier 1 and tier 2)
CGC_v91<- read.csv("cancer_gene_census.csv",row.names = 1)
CGC_genes_v91<- rownames(CGC_v91)
saveRDS(CGC_genes_v91,file = "CGC_v91.rds")

