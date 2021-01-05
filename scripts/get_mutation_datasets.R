# Load mutation data
TCGA_maf<- readRDS("data/TCGA_maf_MC3_selection.rds")

# Turn barcodes into factors (to also retrieve samples without muts later)
TCGA_maf$Tumor_Sample_Barcode<- factor(TCGA_maf$Tumor_Sample_Barcode)

# Add mut_id
TCGA_maf$mut_id<- paste0(TCGA_maf$Hugo_Symbol,"_",gsub("p\\.","",TCGA_maf$HGVSp_Short))

########################################################################
# Select Non-recurrent mutations: 1) Missense, 2) in CGC, 3) NOT Recurrent 
########################################################################

CGC_genes<- readRDS("downloads/cosmic/CGC/v91/CGC_v91.rds")
mut_t<- table(TCGA_maf[TCGA_maf$Variant_Classification=="Missense_Mutation"&TCGA_maf$Hugo_Symbol%in%CGC_genes,"mut_id"])
mut_t<- sort(mut_t[mut_t==1],decreasing=T)
driver_muts<- sample(names(mut_t),500) # 500 random ones

# Select these mutation data
TCGA_maf_driver<- subset(TCGA_maf,TCGA_maf$mut_id%in%driver_muts) 
TCGA_maf_driver$mut_id<- factor(TCGA_maf_driver$mut_id,levels = driver_muts)

# Save
save(TCGA_maf_driver,file="temp/TCGA_maf_noRecurr_driver.RData")

##################################################################################
# Select recurrent no-driver mutations: 1) Missense, 2) NOT in CGC, 3) Recurrent 
#####################################################################################

mut_t<- table(TCGA_maf[TCGA_maf$Variant_Classification=="Missense_Mutation"&!TCGA_maf$Hugo_Symbol%in%CGC_genes,"mut_id"])
mut_t<- sort(mut_t[mut_t>=5],decreasing=T) # 3435
driver_muts<- sample(names(mut_t),500) # 500 random ones

# Select these mutation data
TCGA_maf_driver<- subset(TCGA_maf,TCGA_maf$mut_id%in%driver_muts) 
TCGA_maf_driver$mut_id<- factor(TCGA_maf_driver$mut_id,levels = driver_muts)

# Save
save(TCGA_maf_driver,file="temp/TCGA_maf_recurr_noDriver.RData")

##################################################################################
# Select rnon-ecurrent no-driver mutations: 1) Missense, 2) NOT in CGC, 3) NOT Recurrent 
#####################################################################################

mut_t<- table(TCGA_maf[TCGA_maf$Variant_Classification=="Missense_Mutation"&!TCGA_maf$Hugo_Symbol%in%CGC_genes,"mut_id"])
mut_t<- sort(mut_t[mut_t==1],decreasing=T)
driver_muts<- sample(names(mut_t),500) # 500 random ones

# Select these mutation data
TCGA_maf_driver<- subset(TCGA_maf,TCGA_maf$mut_id%in%driver_muts) 
TCGA_maf_driver$mut_id<- factor(TCGA_maf_driver$mut_id,levels = driver_muts)

# Save
save(TCGA_maf_driver,file="temp/TCGA_maf_noRecurr_noDriver.RData")
