# Download
# wget https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc -O downloads/TCGA_mut/MC3.maf

# Read data
maf_MC3<- read.table("downloads/TCGA_mut/MC3.maf",header = T,stringsAsFactors = F)

# Save data
saveRDS(maf_MC3, "downloads/TCGA_mut/TCGA_maf_MC3.rds")
