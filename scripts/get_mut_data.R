# Load data
maf_MC3 <- readRDS("downloads/TCGA/maf/mc3/TCGA_maf_MC3.rds")

# Select variables we need
TCGA_maf_sel<- maf_MC3[,c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "HGVSp_Short",   "Gene", "CDS_position", "Protein_position", "Amino_acids","ENSP","STRAND")]

# Rename sample barcode
TCGA_maf_sel$Tumor_Sample_Barcode<- substr(TCGA_maf_sel$Tumor_Sample_Barcode,1,15)

# Save data
saveRDS(TCGA_maf_sel,file="data/TCGA_maf_MC3_selection.rds")
