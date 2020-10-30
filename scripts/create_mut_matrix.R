# Load mutation data
TCGA_maf<- readRDS("data/TCGA_maf_MC3_selection.rds")

# Turn barcodes into factors (to also retrieve samples without muts later)
TCGA_maf$Tumor_Sample_Barcode<- factor(TCGA_maf$Tumor_Sample_Barcode)

# Select driver mutations: 1) Missense, 2) in CGC, 3) Recurrent x5
CGC_genes<- readRDS("downloads/cosmic/CGC/v91/CGC_v91.rds")
mut_t<- table(TCGA_maf[TCGA_maf$Variant_Classification=="Missense_Mutation"&TCGA_maf$Hugo_Symbol%in%CGC_genes,"mut_id"])
mut_t<- sort(mut_t[mut_t>=5],decreasing=T) # 688
driver_muts<- names(mut_t)

# Select these mutation data
TCGA_maf_driver<- subset(TCGA_maf,TCGA_maf$mut_id%in%driver_muts) 
TCGA_maf_driver$mut_id<- factor(TCGA_maf_driver$mut_id,levels = driver_muts)

# Create matrix
mut_matrix<- table(TCGA_maf_driver$Tumor_Sample_Barcode,TCGA_maf_driver$mut_id)
# dim(mut_matrix) # 10295x688

# Sort driver muts by frequency
mut_matrix<- mut_matrix[,order(colSums(mut_matrix),decreasing = T)]
# colSums(mut_matrix[,1:5])
# BRAFV600E    IDH1R132H  PIK3CAE545K PIK3CAH1047R     KRASG12D 
# 588          396          277          259          206 
# More or less similar to Marty daya

# Save
saveRDS(mut_matrix,file="data/mut_matrix.rds")

# Also save maf to be used later
save(TCGA_maf_driver,file="temp/TCGA_maf_driver.RData")

