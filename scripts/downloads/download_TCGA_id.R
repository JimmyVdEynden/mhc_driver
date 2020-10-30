# Downloaded together with maf files (mutect2, v10) from GDC data portal at https://portal.gdc.cancer.gov/repository
# After selection and download of manifest, click cases - clinical and download tsv manually

system("tar -xvf clinical.cases_selection.2020-04-03.tar.gz")

# Clinical
TCGA_clin<- read.table("clinical.tsv",header = T,sep="\t")

# Extract cancer ids
TCGA_cancer_id<- as.character(TCGA_clin$project_id)
names(TCGA_cancer_id)<- TCGA_clin$submitter_id
TCGA_cancer_id<- gsub("TCGA-","",TCGA_cancer_id)

# Exposure data: limited smoking, weight and height information
TCGA_exp<- read.table("exposure.tsv",header = T,sep="\t")

# save
saveRDS(TCGA_clin,"TCGA_clin.rds")
saveRDS(TCGA_exp,"TCGA_exp.rds")
saveRDS(TCGA_cancer_id,"TCGA_cancer_id.rds")